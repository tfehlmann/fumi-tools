#include <fumi_tools/dedup.hpp>

#include <algorithm>
#include <random>

#include <nonstd/optional.hpp>
#include <nonstd/string_view.hpp>

#include <cpg/cpg.hpp>

#include <fmt/format.h>

#include <robin_hood/robin_hood.h>

#include <htslib/sam.h>

#include <fumi_tools/cast_helper.hpp>
#include <fumi_tools/helper.hpp>
#include <fumi_tools/umi_clusterer.hpp>

namespace fumi_tools {

namespace {
std::mt19937 rand_gen;
std::uniform_real_distribution<> udistrib(0, 1);
}  // namespace

nonstd::string_view get_umi(const char* c, uint8_t len) {
  auto res = nonstd::string_view(c, len);
  auto pos = res.rfind('_');
  if (pos == nonstd::string_view::npos) {
    throw std::runtime_error(fmt::format("Did not find umi for read {}!", res.to_string()));
  }
  return res.substr(res.rfind('_') + 1);
}

/**
 * Takes a cigar string and finds the first splice position as
    an offset from the start. To find the 5' end (read coords) of
    the junction for a reverse read, pass in the reversed cigar tuple
 */
uint32_t find_splice(uint32_t* cigar, unsigned int n_cigar, bool reverse) {
  uint32_t offset = 0;

  auto cigar_end = !reverse ? as_signed(n_cigar) : -1l;
  auto cigar_start = !reverse ? 0l : as_signed(n_cigar) - 1l;
  // a soft clip at the end of the read is taken as splicing
  // where as a soft clip at the start is not.
  if ((cigar[cigar_start] & BAM_CIGAR_MASK) == BAM_CSOFT_CLIP) {
    offset = cigar[cigar_start] & BAM_CIGAR_SHIFT;
    if (reverse) {
      --cigar_start;
    } else {
      ++cigar_start;
    }
  }

  for (; cigar_start != cigar_end; reverse ? --cigar_start : ++cigar_start) {
    switch (cigar[cigar_start] & BAM_CIGAR_MASK) {
      case BAM_CSOFT_CLIP:
      case BAM_CREF_SKIP:
        return offset;
      case BAM_CMATCH:
      case BAM_CDEL:
      case BAM_CEQUAL:
      case BAM_CDIFF:
        offset += cigar[cigar_start] & BAM_CIGAR_SHIFT;
        break;
    }
  }
  return 0;
}

bool cigar_has_cref_skip(uint32_t* cigar, unsigned int len) {
  for (auto i = 0ul; i < len; ++i) {
    if ((cigar[i] & BAM_CIGAR_MASK) == BAM_CREF_SKIP) {
      return true;
    }
  }
  return false;
}

std::tuple<uint64_t, uint64_t, bool> get_read_position(
    const bam1_t* read, uint32_t soft_clip_threshold) {
  auto is_spliced = false;
  auto* cigar = bam_get_cigar(read);
  auto n_cigar = read->core.n_cigar;
  if ((read->core.flag & BAM_FUNMAP) != 0) {
    auto pos = bam_endpos(read);
    if ((cigar[n_cigar - 1] & BAM_CIGAR_MASK) == BAM_CSOFT_CLIP) {
      auto count = cigar[n_cigar - 1] >> BAM_CIGAR_SHIFT;
      pos += count;
    }
    auto start = read->core.pos;
    if (cigar_has_cref_skip(cigar, n_cigar) ||
        ((cigar[0] & BAM_CIGAR_MASK) == BAM_CSOFT_CLIP &&
         (cigar[0] & BAM_CIGAR_SHIFT) > soft_clip_threshold)) {
      is_spliced = find_splice(cigar, n_cigar, true);
    }
    return std::tie(start, pos, is_spliced);
  } else {
    auto pos = read->core.pos;
    if ((cigar[0] & BAM_CIGAR_MASK) == BAM_CSOFT_CLIP) {
      pos -= cigar[0] & BAM_CIGAR_SHIFT;
    }
    auto start = pos;
    if (cigar_has_cref_skip(cigar, n_cigar) ||
        ((cigar[n_cigar - 1] & BAM_CIGAR_MASK) == BAM_CSOFT_CLIP &&
         (cigar[n_cigar - 1] & BAM_CIGAR_SHIFT) > soft_clip_threshold)) {
      is_spliced = find_splice(cigar, n_cigar, false);
    }
    return std::tie(start, pos, is_spliced);
  }
}

void update_read_map(
    const bam1_t* read, uint64_t pos, read_group key, std::string umi,
    robin_hood::unordered_flat_map<
        uint64_t,
        robin_hood::unordered_flat_map<
            read_group,
            robin_hood::unordered_flat_map<
                std::string,
                std::pair<std::unique_ptr<bam1_t, bam1_t_deleter>, uint64_t>>>>&
        read_map,
    robin_hood::unordered_flat_map<
        uint64_t,
        robin_hood::unordered_flat_map<
            read_group, robin_hood::unordered_flat_map<std::string, uint64_t>>>&
        read_counts) {
  auto it = read_map[pos][key].find(umi);
  if (it != std::end(read_map[pos][key])) {
    auto& res = it->second;
    res.second += 1;
    auto read_qual = bam_get_qual(read);
    auto other_qual = bam_get_qual(res.first);
    if (read_qual < other_qual) {
      return;
    }
    if (read_qual > other_qual) {
      res.first.reset(bam_dup1(read));
      read_counts[pos][key][umi] = 0;
    }
    auto& count_res = read_counts[pos][key][umi];
    ++count_res;
    auto prob = 1.0 / count_res;

    if (udistrib(rand_gen) < prob) {
      res.first.reset(bam_dup1(read));
    }
  } else {
    read_map[pos][key][umi] = {
        std::unique_ptr<bam1_t, bam1_t_deleter>(bam_dup1(read)), 1ul};
    read_counts[pos][key][umi] = 0;
  }
}

template <class Fun>
void process_bam_read_chunks(samFile* file, bam_hdr_t* bam_hdr, umi_opts opts,
                             Fun fun) {
  auto cur_ref = 0;
  auto last_ref = -1;
  auto last_pos = 0ul;
  auto last_output_pos = 0ul;

  robin_hood::unordered_flat_map<
      uint64_t,
      robin_hood::unordered_flat_map<
          read_group,
          robin_hood::unordered_flat_map<
              std::string,
              std::pair<std::unique_ptr<bam1_t, bam1_t_deleter>, uint64_t>>>>
      read_map;
  robin_hood::unordered_flat_map<
      uint64_t,
      robin_hood::unordered_flat_map<
          read_group, robin_hood::unordered_flat_map<std::string, uint64_t>>>
      read_counts;

  auto output_positions = [&read_map, &read_counts,
                           &fun](nonstd::optional<uint64_t> start) {
    std::vector<uint64_t> positions;
    positions.reserve(read_map.size());
    for (auto& k_v : read_map) {
      if (!start.has_value() || k_v.first + 1000 < start) {
        positions.push_back(k_v.first);
      }
    }
    std::sort(positions.begin(), positions.end());
    for (auto p : positions) {
      std::vector<read_group> sorted_keys;
      sorted_keys.reserve(read_map.size());
      auto& map_p = read_map[p];
      for (auto& k_v : read_map[p]) {
        sorted_keys.push_back(k_v.first);
      }
      std::sort(sorted_keys.begin(), sorted_keys.end());
      for (auto& k : sorted_keys) {
        fun(map_p[k]);
      }
    }
    for (auto p : positions) {
      read_map.erase(p);
      read_counts.erase(p);
    }
  };

  cpg::cpg_cfg prog_cfg{};
  prog_cfg.unit = "aln";
  prog_cfg.unit_scale = true;
  prog_cfg.mininterval = 3;

  auto progress = cpg::cpg(prog_cfg);

  bam1_t* record = bam_init1();
  while (sam_read1(file, bam_hdr, record) > 0) {
    if ((record->core.flag & BAM_FUNMAP) == 0) {
      cur_ref = record->core.tid;
      auto* qname = bam_get_qname(record);
      auto umi = get_umi(qname, ::strlen(qname));
      uint64_t start = 0;
      uint64_t pos = 0;
      bool is_spliced = false;
      std::tie(start, pos, is_spliced) =
          get_read_position(record, opts.soft_clip_threshold);

      // new ref, so output all previous reads
      if(cur_ref != last_ref){
          output_positions(nonstd::nullopt);
          last_output_pos = 0;
      } else if (last_output_pos + 1000 < start) {
        output_positions(start);
        last_output_pos = start;
      }

      last_pos = std::max(pos, start);
      last_ref = cur_ref;
      auto key = read_group{
          bam_is_rev(record), opts.spliced && is_spliced != 0,
          static_cast<uint16_t>(opts.read_length ? record->core.l_qseq : 0)};
      update_read_map(record, pos, key, std::string{umi}, read_map,
                      read_counts);
    }
    progress.update();
  }
  output_positions(nonstd::nullopt);
  bam_destroy1(record);
}

void dedup(const std::string& input, const std::string& output, umi_opts opts) {
  rand_gen.seed(opts.seed);

  samFile* file = hts_open(input.c_str(), "r");

  if (file == nullptr) {
    throw std::runtime_error(fmt::format("Could not open file '{}'", input));
  }

  // read header
  bam_hdr_t* bam_hdr = sam_hdr_read(file);

  // require sorted bam file
  if ((bam_hdr->l_text > 3) && (strncmp(bam_hdr->text, "@HD", 3) == 0)) {
    auto* p = strstr(bam_hdr->text, "\tSO:coordinate");
    auto* q = strchr(bam_hdr->text, '\n');
    // Looking for SO:coordinate within the @HD line only
    // (e.g. must ignore in a @CO comment line later in header)
    if ((p == nullptr) || (p >= q)) {
      throw std::runtime_error(
          fmt::format("BAM file needs to be coordinate sorted!"));
    }
  } else {
    throw std::runtime_error(
        fmt::format("BAM file needs to be coordinate sorted!"));
  }

  samFile* out = hts_open(output.c_str(), "wb");
  if (out == nullptr) {
    throw std::runtime_error(fmt::format("Could not open file '{}'", output));
  }
  if (sam_hdr_write(out, bam_hdr) < -1) {
    throw std::runtime_error(
        fmt::format("Could not write header to file '{}'", output));
  }

  umi_clusterer clusterer(opts.method);
  process_bam_read_chunks(
      file, bam_hdr, opts, [&clusterer, &out, &bam_hdr](auto& bundle) {
        clusterer(bundle, [&out, &bam_hdr](auto& read, auto& umi, auto& count) {
          sam_write1(out, bam_hdr, read.get());
        });
      });
  bam_hdr_destroy(bam_hdr);

  hts_close(file);
  hts_close(out);
}

}  // namespace fumi_tools
