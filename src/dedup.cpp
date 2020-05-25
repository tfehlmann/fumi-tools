#include <fmt/format.h>
#include <htslib/sam.h>
#include <robin_hood/robin_hood.h>

#include <algorithm>
#include <cpg/cpg.hpp>
#include <fumi_tools/cast_helper.hpp>
#include <fumi_tools/dedup.hpp>
#include <fumi_tools/helper.hpp>
#include <fumi_tools/umi_clusterer.hpp>
#include <iostream>
#include <nonstd/optional.hpp>
#include <nonstd/string_view.hpp>
#include <random>

#define SHOW_DEBUG_OUTPUT false

namespace fumi_tools {

namespace {
std::mt19937 rand_gen;
std::uniform_real_distribution<> udistrib(0, 1);
}  // namespace

nonstd::string_view get_umi(const char* c, uint8_t len) {
  auto res = nonstd::string_view(c, len);
  auto pos = res.rfind('_');
  if (pos == nonstd::string_view::npos) {
    throw std::runtime_error(
        fmt::format("Did not find umi for read {}!", res.to_string()));
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
    offset = cigar[cigar_start] >> BAM_CIGAR_SHIFT;
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
        offset += cigar[cigar_start] >> BAM_CIGAR_SHIFT;
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
    const bam1_t* read,
    uint32_t soft_clip_threshold) {
  auto is_spliced = false;
  auto* cigar = bam_get_cigar(read);
  auto n_cigar = read->core.n_cigar;
  if ((read->core.flag & BAM_FREVERSE) != 0) {
    auto pos = bam_endpos(read);
    if ((cigar[n_cigar - 1] & BAM_CIGAR_MASK) == BAM_CSOFT_CLIP) {
      auto count = cigar[n_cigar - 1] >> BAM_CIGAR_SHIFT;
      pos += count;
    }
    auto start = read->core.pos;
    if (cigar_has_cref_skip(cigar, n_cigar) ||
        ((cigar[0] & BAM_CIGAR_MASK) == BAM_CSOFT_CLIP &&
         (cigar[0] >> BAM_CIGAR_SHIFT) > soft_clip_threshold)) {
      is_spliced = find_splice(cigar, n_cigar, true);
    }
    return std::make_tuple(start, pos, is_spliced);
  } else {
    auto pos = read->core.pos;
    if ((cigar[0] & BAM_CIGAR_MASK) == BAM_CSOFT_CLIP) {
      pos -= cigar[0] >> BAM_CIGAR_SHIFT;
    }
    auto start = pos;
    if (cigar_has_cref_skip(cigar, n_cigar) ||
        ((cigar[n_cigar - 1] & BAM_CIGAR_MASK) == BAM_CSOFT_CLIP &&
         (cigar[n_cigar - 1] >> BAM_CIGAR_SHIFT) > soft_clip_threshold)) {
      is_spliced = find_splice(cigar, n_cigar, false);
    }
    return std::make_tuple(start, pos, is_spliced);
  }
}

struct custom_bam1_hash {
  template <typename deleter>
  std::size_t operator()(const std::unique_ptr<bam1_t, deleter>& lhs) const {
    return this->operator()(lhs.get());
  }

  std::size_t operator()(const bam1_t* lhs) const {
    auto hi_lhs = bam_aux_get(lhs, "HI");
    auto hi_lhs_int = hi_lhs == nullptr ? 0 : bam_aux2i(hi_lhs);
    return robin_hood::hash_int(lhs->core.pos) ^
           robin_hood::hash_int(lhs->core.tid) ^
           robin_hood::hash_bytes(bam_get_qname(lhs),
                                  std::strlen(bam_get_qname(lhs))) ^
           robin_hood::hash_int(lhs->core.mpos) ^
           robin_hood::hash_int(lhs->core.mtid) ^
           robin_hood::hash_int(lhs->core.isize) ^
           robin_hood::hash_int(hi_lhs_int);
  }
};

struct custom_bam1_eq {
  template <typename deleter>
  bool operator()(const std::unique_ptr<bam1_t, deleter>& lhs,
                  const std::unique_ptr<bam1_t, deleter>& rhs) const {
    return this->operator()(lhs.get(), rhs.get());
  }

  bool operator()(const bam1_t* lhs, const bam1_t* rhs) const {
    auto check1 =
        lhs->core.pos == rhs->core.pos && lhs->core.tid == rhs->core.tid &&
        lhs->core.mpos == rhs->core.mpos && lhs->core.mtid == rhs->core.mtid &&
        lhs->core.isize == rhs->core.isize;
    if (!check1 || nonstd::string_view(bam_get_qname(lhs),
                                       std::strlen(bam_get_qname(lhs))) !=
                       nonstd::string_view(bam_get_qname(rhs),
                                           std::strlen(bam_get_qname(rhs)))) {
      return false;
    } else {
      auto hi_lhs = bam_aux_get(lhs, "HI");
      bool same_hi = true;
      if (hi_lhs != nullptr) {
        auto hi_rhs = bam_aux_get(rhs, "HI");
        if (hi_rhs) {
          same_hi = bam_aux2i(hi_lhs) == bam_aux2i(hi_rhs);
        }
      }
      return same_hi;
    }
  }
};

bam1_t build_mate_bam1_dummy(const bam1_t& read) {
  bam1_t dummy;
  dummy.core.pos = read.core.mpos;
  dummy.core.tid = read.core.mtid;
  dummy.core.n_cigar = read.core.n_cigar;
  dummy.core.l_qname = read.core.l_qname;
  dummy.core.l_qseq = read.core.l_qseq;
  dummy.core.mpos = read.core.pos;
  dummy.core.mtid = read.core.tid;
  dummy.core.isize = -read.core.isize;
  dummy.data = read.data;
  dummy.l_data = read.l_data;
  return dummy;
}

bool read_is_potentially_after_mate(const bam1_t& read){
  return (read.core.mtid < read.core.tid ||
   (read.core.mtid == read.core.tid && read.core.pos >= read.core.mpos));
}

bool read_is_after_mate(const bam1_t& read){
  return (read.core.mtid < read.core.tid ||
   (read.core.mtid == read.core.tid && read.core.pos > read.core.mpos));
}

template <class ReadGroup, bool is_paired>
void update_read_map(
    bam1_t* read,
    uint64_t pos,
    ReadGroup key,
    std::string umi,
    robin_hood::unordered_flat_map<
        uint64_t,
        robin_hood::unordered_flat_map<
            ReadGroup,
            robin_hood::unordered_flat_map<
                std::string,
                std::pair<std::unique_ptr<bam1_t, bam1_t_deleter>, uint64_t>>>>&
        read_map,
    robin_hood::unordered_flat_map<
        uint64_t,
        robin_hood::unordered_flat_map<
            ReadGroup,
            robin_hood::unordered_flat_map<std::string, uint64_t>>>&
        read_counts,
    robin_hood::unordered_flat_set<std::unique_ptr<bam1_t, bam1_t_deleter>,
                                   custom_bam1_hash,
                                   custom_bam1_eq>& paired_read_map,
    robin_hood::unordered_set<bam1_t*, custom_bam1_hash, custom_bam1_eq>&
        current_reads) {
  auto it = read_map[pos][key].find(umi);
  if (it != std::end(read_map[pos][key])) {
    auto& res = it->second;
    res.second += 1;
    auto read_qual = read->core.qual;
    auto other_qual = res.first->core.qual;
    if (read_qual < other_qual) {
      // bad qual so drop pair
      if (is_paired && read_is_potentially_after_mate(*read)) {
        bam1_t dummy = build_mate_bam1_dummy(*read);
        std::unique_ptr<bam1_t, bam1_t_deleter> dummy_ptr(&dummy);
        paired_read_map.erase(dummy_ptr);
        dummy_ptr.release();
      }
      return;
    }
    if (read_qual > other_qual) {
      // replace with other read, so remove paired
      if (is_paired && read_is_potentially_after_mate(*res.first)) {
        bam1_t dummy = build_mate_bam1_dummy(*res.first);
        std::unique_ptr<bam1_t, bam1_t_deleter> dummy_ptr(&dummy);
        paired_read_map.erase(dummy_ptr);
        dummy_ptr.release();
      }
      if (is_paired) {
        current_reads.erase(res.first.get());
      }
      auto* new_read = bam_dup1(read);
      res.first.reset(new_read);
      if (is_paired) {
        current_reads.insert(new_read);
      }
      read_counts[pos][key][umi] = 0;
      return;
    }
    auto& count_res = read_counts[pos][key][umi];
    ++count_res;
    auto prob = 1.0 / count_res;

    if (udistrib(rand_gen) < prob) {
      if (is_paired && read_is_potentially_after_mate(*res.first)) {
        bam1_t dummy = build_mate_bam1_dummy(*res.first);
        std::unique_ptr<bam1_t, bam1_t_deleter> dummy_ptr(&dummy);
        paired_read_map.erase(dummy_ptr);
        dummy_ptr.release();
      }
      if (is_paired) {
        current_reads.erase(res.first.get());
      }
      auto* new_read = bam_dup1(read);
      res.first.reset(new_read);
      if (is_paired) {
        current_reads.insert(new_read);
      }
    } else {
      // bad qual so drop pair
      if (is_paired && read_is_potentially_after_mate(*read)) {
        bam1_t dummy = build_mate_bam1_dummy(*read);
        std::unique_ptr<bam1_t, bam1_t_deleter> dummy_ptr(&dummy);
        paired_read_map.erase(dummy_ptr);
        dummy_ptr.release();
      }
    }
  } else {
    auto* new_read = bam_dup1(read);
    if (is_paired) {
      current_reads.insert(new_read);
    }
    read_map[pos][key][umi] = {
        std::unique_ptr<bam1_t, bam1_t_deleter>(new_read), 1ul};
    read_counts[pos][key][umi] = 0;
  }
}

template <class ReadGroup>
std::ostream& operator<<(
    std::ostream& out,
    const robin_hood::unordered_flat_map<
        ReadGroup,
        robin_hood::unordered_flat_map<
            std::string,
            std::pair<std::unique_ptr<bam1_t, bam1_t_deleter>, uint64_t>>>&
        lhs) {
  out << '{';
  for (auto& e : lhs) {
    out << e.first << ':';
    for (auto& el : e.second) {
      out << el.first << ' ';
    }
    out << '\n';
  }
  out << '}';
  return out;
}

std::ostream& operator<<(std::ostream& out,
                         const robin_hood::unordered_flat_set<
                             std::unique_ptr<bam1_t, bam1_t_deleter>,
                             custom_bam1_hash,
                             custom_bam1_eq>& lhs) {
  out << '{';
  for (auto& e : lhs) {
    out << bam_get_qname(e) << '|' << e->core.pos << ",\n";
  }
  out << '}';
  return out;
}

template <class ReadGroup, bool is_paired, class Fun>
void process_bam_read_chunks_helper(samFile* file,
                                    bam_hdr_t* bam_hdr,
                                    umi_opts opts,
                                    samFile* out,
                                    Fun fun) {
  auto cur_ref = 0;
  auto last_ref = -1;
  auto last_pos = 0ul;
  auto last_output_pos = 0ul;

  robin_hood::unordered_flat_map<
      uint64_t,
      robin_hood::unordered_flat_map<
          ReadGroup,
          robin_hood::unordered_flat_map<
              std::string,
              std::pair<std::unique_ptr<bam1_t, bam1_t_deleter>, uint64_t>>>>
      read_map;
  robin_hood::unordered_flat_map<
      uint64_t,
      robin_hood::unordered_flat_map<
          ReadGroup, robin_hood::unordered_flat_map<std::string, uint64_t>>>
      read_counts;

  robin_hood::unordered_set<bam1_t*, custom_bam1_hash, custom_bam1_eq>
      current_reads;

  robin_hood::unordered_flat_set<std::unique_ptr<bam1_t, bam1_t_deleter>,
                                 custom_bam1_hash, custom_bam1_eq>
      not_yet_paired_reads;
  robin_hood::unordered_flat_set<std::unique_ptr<bam1_t, bam1_t_deleter>,
                                 custom_bam1_hash, custom_bam1_eq>
      paired_read_map;

  auto output_positions = [&read_map, &read_counts, &fun, &current_reads,
                           &paired_read_map, &not_yet_paired_reads](
                              nonstd::optional<uint64_t> start,
                              int32_t bam_pos) {
    std::vector<uint64_t> positions;
    positions.reserve(read_map.size());
    for (auto& k_v : read_map) {
      if (!start.has_value() || k_v.first + 1000 < start) {
        positions.push_back(k_v.first);
      }
    }
    std::sort(positions.begin(), positions.end());
    for (auto p : positions) {
      auto& map_p = read_map[p];
      std::vector<ReadGroup> sorted_keys;
      sorted_keys.reserve(map_p.size());
      for (auto& k_v : read_map[p]) {
        sorted_keys.push_back(k_v.first);
      }
      std::sort(sorted_keys.begin(), sorted_keys.end());
      for (auto& k : sorted_keys) {
        if (is_paired) {
          for (auto& e : map_p[k]) {
            current_reads.erase(e.second.first.get());
          }
        }
        fun(map_p[k], bam_pos, paired_read_map, not_yet_paired_reads);
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
  prog_cfg.desc = "Dedup UMIs";

  auto progress = cpg::cpg(prog_cfg);

  bam1_t* record = bam_init1();
  while (sam_read1(file, bam_hdr, record) > 0) {
    if ((record->core.flag & BAM_FUNMAP) == 0) {
      if (is_paired && (record->core.flag & BAM_FMUNMAP) != 0 &&
          opts.unpaired_reads == "discard") {
        continue;
      }
      if (is_paired && (record->core.flag & BAM_FREAD2) != 0) {
        if (record->core.tid <= cur_ref) {
          // we already saw r1
          if (record->core.mpos < record->core.pos || record->core.tid < cur_ref) {
            bam1_t dummy = build_mate_bam1_dummy(*record);
            // keep paired read only if we kept r1
            if (current_reads.find(&dummy) != current_reads.end()) {
              paired_read_map.emplace(bam_dup1(record));
            } else {
              // maybe r1 is from a previous bundle
              std::unique_ptr<bam1_t, bam1_t_deleter> dummy_ptr(&dummy);
              auto it = not_yet_paired_reads.find(dummy_ptr);
              if (it != not_yet_paired_reads.end()) {
                // output directly
                if (sam_write1(out, bam_hdr, it->get()) < 0) {
                  std::cerr << "Failed to write to output file!" << std::endl;
                  std::exit(1);
                }
                if (sam_write1(out, bam_hdr, record) < 0) {
                  std::cerr << "Failed to write to output file!" << std::endl;
                  std::exit(1);
                }
                not_yet_paired_reads.erase(it);
              }
              // this is not a pointer to the free store, so release ownership
              dummy_ptr.release();
              continue;
            }
          } else {  // >= 0, need to check later
            paired_read_map.emplace(bam_dup1(record));
          }
        } else {
          paired_read_map.emplace(bam_dup1(record));
        }
        continue;
      }
      cur_ref = record->core.tid;
      auto* qname = bam_get_qname(record);
      auto umi = get_umi(qname, std::strlen(qname));
      uint64_t start = 0;
      uint64_t pos = 0;
      bool is_spliced = false;
      std::tie(start, pos, is_spliced) =
          get_read_position(record, opts.soft_clip_threshold);

      if (is_paired && (record->core.flag & BAM_FPAIRED) != 0 &&
          record->core.tid != record->core.mtid) {
        // chimeric read pair
        // other possibility is "use", which is implicitly handled
        if (opts.chimeric_pairs == "discard") {
          bam1_t dummy = build_mate_bam1_dummy(*record);
          std::unique_ptr<bam1_t, bam1_t_deleter> dummy_ptr(&dummy);
          paired_read_map.erase(dummy_ptr);
          dummy_ptr.release();
          continue;
        }
      }

      // new ref, so output all previous reads
      if (cur_ref != last_ref) {
        output_positions(nonstd::nullopt, std::numeric_limits<int32_t>::max());
        last_output_pos = 0;
      } else if (last_output_pos + 1000 < start) {
        output_positions(start, record->core.pos);
        last_output_pos = start;
      }

      last_pos = std::max(pos, start);
      last_ref = cur_ref;
      auto key = ReadGroup(
          bam_is_rev(record), opts.spliced && is_spliced != 0,
          (!opts.ignore_tlen && opts.paired) ? record->core.isize : 0,
          static_cast<uint16_t>(opts.read_length ? record->core.l_qseq : 0));
      update_read_map<ReadGroup, is_paired>(record, pos, key, std::string{umi},
                                            read_map, read_counts,
                                            paired_read_map, current_reads);
    }
    progress.update();
  }
  output_positions(nonstd::nullopt, std::numeric_limits<int32_t>::max());
  if (is_paired && SHOW_DEBUG_OUTPUT && !not_yet_paired_reads.empty()) {
    std::cerr << not_yet_paired_reads << std::endl;
  }
  if (is_paired && opts.unpaired_reads == "use") {
    for (auto& read : not_yet_paired_reads) {
      if (sam_write1(out, bam_hdr, read.get()) < 0) {
        std::cerr << "Failed to write to output file!" << std::endl;
        std::exit(1);
      }
    }
  }
  bam_destroy1(record);
}

template <class Fun>
void process_bam_read_chunks(samFile* file,
                             bam_hdr_t* bam_hdr,
                             umi_opts opts,
                             samFile* out,
                             Fun fun) {
  if (opts.paired) {
    process_bam_read_chunks_helper<read_group_paired, true>(file, bam_hdr, opts,
                                                            out, fun);
  } else {
    process_bam_read_chunks_helper<read_group, false>(file, bam_hdr, opts, out,
                                                      fun);
  }
}

template <class Fun>
void process_paired_reads(
    std::unique_ptr<bam1_t, bam1_t_deleter>& r1,
    robin_hood::unordered_flat_set<std::unique_ptr<bam1_t, bam1_t_deleter>,
                                   custom_bam1_hash,
                                   custom_bam1_eq>& paired_reads,
    robin_hood::unordered_flat_set<std::unique_ptr<bam1_t, bam1_t_deleter>,
                                   custom_bam1_hash,
                                   custom_bam1_eq>& not_yet_paired_reads,
    int32_t bam_pos,
    Fun fun) {
  if ((r1->core.flag & BAM_FMUNMAP) == 0) {  // has mate
    if (r1->core.mpos <= bam_pos) {  // we can only have the mate if it was
                                     // before our current position
      bam1_t dummy = build_mate_bam1_dummy(*r1);
      std::unique_ptr<bam1_t, bam1_t_deleter> dummy_ptr(&dummy);
      auto it = paired_reads.find(dummy_ptr);
      dummy_ptr.release();
      if (it != paired_reads.end()) {
        fun(r1.get());
        fun(it->get());
        paired_reads.erase(it);
      } else {
        not_yet_paired_reads.insert(std::move(r1));
      }
    } else {
      not_yet_paired_reads.insert(std::move(r1));
    }
  } else {
    fun(r1.get());
  }
}

void dedup(const std::string& input, const std::string& output, umi_opts opts) {
  rand_gen.seed(opts.seed);

  samFile* file = hts_open(input.c_str(), "r");

  if (file == nullptr) {
    throw std::runtime_error(fmt::format("Could not open file '{}'", input));
  }

  if (opts.ithreads > 1) {
    hts_set_threads(file, opts.ithreads);
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

  samFile* out = hts_open(
      output.c_str(),
      opts.uncompressed ? "wbu" : ends_with(output, ".bam") ? "wb" : "w");
  if (out == nullptr) {
    throw std::runtime_error(fmt::format("Could not open file '{}'", output));
  }

  if (opts.othreads > 1) {
    hts_set_threads(out, opts.othreads);
  }

  if (sam_hdr_write(out, bam_hdr) < -1) {
    throw std::runtime_error(
        fmt::format("Could not write header to file '{}'", output));
  }

  umi_clusterer clusterer(opts.method);

  process_bam_read_chunks(
      file, bam_hdr, opts, out,
      [&clusterer, &out, &bam_hdr, paired = opts.paired](
          auto& bundle, auto bam_pos, auto& paired_read_map,
          auto& not_yet_paired_reads) {
        clusterer(bundle, [&out, &bam_hdr, &paired_read_map,
                           &not_yet_paired_reads, &bam_pos,
                           paired](auto& read, auto& /*umi*/, auto& /*count*/) {
          if (paired) {
            process_paired_reads(
                read, paired_read_map, not_yet_paired_reads, bam_pos,
                [&out, &bam_hdr](const bam1_t* read) {
                  if (sam_write1(out, bam_hdr, read) < 0) {
                    std::cerr << "Failed to write to output file!" << std::endl;
                    std::exit(1);
                  }
                });
          } else {
            if (sam_write1(out, bam_hdr, read.get()) < 0) {
              std::cerr << "Failed to write to output file!" << std::endl;
              std::exit(1);
            }
          }
        });
      });
  bam_hdr_destroy(bam_hdr);

  hts_close(file);
  hts_close(out);
}

}  // namespace fumi_tools
