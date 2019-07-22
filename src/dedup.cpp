#include <fumi_tools/dedup.hpp>

#include <algorithm>
#include <random>
#include <string_view>

#include <fmt/format.h>

#include <robin_hood/robin_hood.h>

#include <htslib/sam.h>


#include <seqan/bam_io.h>
#include <seqan/parallel.h>
#include <seqan/parallel/parallel_macros.h>



#include <fumi_tools/cast_helper.hpp>
#include <fumi_tools/umi_clusterer.hpp>
#include <fumi_tools/helper.hpp>


namespace fumi_tools {


namespace {
std::mt19937 rand_gen;
std::uniform_real_distribution<> udistrib(0, 1);
}  // namespace

std::string_view get_umi(const seqan::CharString& c) {
  auto res = std::string_view(seqan::toCString(c), seqan::length(c));
  auto pos = res.rfind('_');
  if (pos == std::string_view::npos) {
    throw std::runtime_error(fmt::format("Did not find umi for read {}!", res));
  }
  return res.substr(res.rfind('_') + 1);
}

std::string_view get_umi(const char* c, uint8_t len) {
  auto res = std::string_view(c, len);
  auto pos = res.rfind('_');
  if (pos == std::string_view::npos) {
    throw std::runtime_error(fmt::format("Did not find umi for read {}!", res));
  }
  return res.substr(res.rfind('_') + 1);
}

/**
 * Takes a cigar string and finds the first splice position as
    an offset from the start. To find the 5' end (read coords) of
    the junction for a reverse read, pass in the reversed cigar tuple
 */
template <class CigarString>
uint32_t find_splice(CigarString& cigar, bool reverse) {
  uint32_t offset = 0;

  auto cigar_end = !reverse ? as_signed(seqan::length(cigar)) : -1l;
  auto cigar_start = !reverse ? 0l : as_signed(seqan::length(cigar)) - 1l;
  // a soft clip at the end of the read is taken as splicing
  // where as a soft clip at the start is not.
  if (cigar[cigar_start].operation == 'S') {
    offset = cigar[cigar_start].count;
    if (reverse) {
      --cigar_start;
    } else {
      ++cigar_start;
    }
  }

  for (; cigar_start != cigar_end; reverse ? --cigar_start : ++cigar_start) {
    switch (cigar[cigar_start].operation) {
      case 'N':
      case 'S':
        return offset;
      case 'M':
      case 'D':
      case '=':
      case 'X':
        offset += cigar[cigar_start].count;
        break;
      case 'I':
      case 'H':
      case 'P':
        continue;
    }
  }
  return 0;
}

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

std::tuple<uint64_t, uint64_t, bool> get_read_position(
    const seqan::BamAlignmentRecord& read, uint32_t soft_clip_threshold) {
  auto is_spliced = false;
  if (seqan::hasFlagRC(read)) {
    auto aln_len = seqan::getAlignmentLengthInRef(read);
    auto pos = read.beginPos + aln_len;
    if (read.cigar[seqan::length(read.cigar) - 1].operation == 'S') {
      pos += read.cigar[seqan::length(read.cigar) - 1].count;
    }
    auto start = read.beginPos;
    if (std::any_of(seqan::begin(read.cigar), seqan::end(read.cigar),
                    [](auto& e) { return e.operation == 'N'; }) ||
        (read.cigar[0].operation == 'S' &&
         read.cigar[0].count > soft_clip_threshold)) {
      is_spliced = find_splice(read.cigar, true);
    }
    return std::tie(start, pos, is_spliced);
  } else {
    auto pos = read.beginPos;
    if (read.cigar[0].operation == 'S') {
      pos -= read.cigar[0].count;
    }
    auto start = pos;
    if (std::any_of(seqan::begin(read.cigar), seqan::end(read.cigar),
                    [](auto& e) { return e.operation == 'N'; }) ||
        (read.cigar[seqan::length(read.cigar) - 1].operation == 'S' &&
         read.cigar[seqan::length(read.cigar) - 1].count >
             soft_clip_threshold)) {
      is_spliced = find_splice(read.cigar, false);
    }
    return std::tie(start, pos, is_spliced);
  }
}

bool cigar_has_cref_skip(uint32_t* cigar, unsigned int len){
    for(auto i = 0ul; i < len; ++i){
        if((cigar[i] & BAM_CIGAR_MASK) == BAM_CREF_SKIP){
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
    //auto aln_len = seqan::getAlignmentLengthInRef(read);
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
    if ((cigar[0] & BAM_CIGAR_MASK) == 'S') {
      pos -= cigar[0] & BAM_CIGAR_SHIFT;
    }
    auto start = pos;
    if (cigar_has_cref_skip(cigar, n_cigar) ||
        ((cigar[n_cigar - 1] & BAM_CIGAR_MASK) == BAM_CSOFT_CLIP &&
         (cigar[n_cigar - 1] & BAM_CIGAR_SHIFT) >
             soft_clip_threshold)) {
      is_spliced = find_splice(cigar, n_cigar, false);
    }
    return std::tie(start, pos, is_spliced);
  }
}


void update_read_map(
    const seqan::BamAlignmentRecord& read, uint64_t pos,
    read_group key, std::string umi,
    robin_hood::unordered_flat_map<
        uint64_t,
        robin_hood::unordered_flat_map<
            read_group,
            robin_hood::unordered_flat_map<
                std::string, std::pair<seqan::BamAlignmentRecord, uint64_t>>>>&
        read_map,
    robin_hood::unordered_flat_map<
        uint64_t, robin_hood::unordered_flat_map<
                      read_group,
                      robin_hood::unordered_flat_map<std::string, uint64_t>>>&
        read_counts) {
  auto it = read_map[pos][key].find(umi);
  if (it != std::end(read_map[pos][key])) {
    auto& res = it->second;
    res.second += 1;
    if (read.mapQ < res.first.mapQ) {
      return;
    }
    if (read.mapQ > res.first.mapQ) {
      res.first = read;
      read_counts[pos][key][umi] = 0;
    }
    auto& count_res = read_counts[pos][key][umi];
    ++count_res;
    auto prob = 1.0 / count_res;

    if (udistrib(rand_gen) < prob) {
      res.first = read;
    }
  } else {
    read_map[pos][key][umi] = {read, 1ul};
    read_counts[pos][key][umi] = 0;
  }
}

void update_read_map(
    const bam1_t* read, uint64_t pos,
    read_group key, std::string umi,
    robin_hood::unordered_flat_map<
        uint64_t,
        robin_hood::unordered_flat_map<
            read_group,
            robin_hood::unordered_flat_map<
                std::string, std::pair<std::unique_ptr<bam1_t, bam1_t_deleter>, uint64_t>>>>&
        read_map,
    robin_hood::unordered_flat_map<
        uint64_t, robin_hood::unordered_flat_map<
                      read_group,
                      robin_hood::unordered_flat_map<std::string, uint64_t>>>&
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
    read_map[pos][key][umi] = {std::unique_ptr<bam1_t, bam1_t_deleter>(bam_dup1(read)), 1ul};
    read_counts[pos][key][umi] = 0;
  }
}

template<class Fun>
void process_bam_read_chunks(/*seqan::BamFileIn&*/ samFile* file, bam_hdr_t* bam_hdr, umi_opts opts, Fun fun){
    auto cur_ref = 0;
    auto last_ref = -1;
    auto last_pos = 0ul;
    auto last_output_pos = 0ul;

    robin_hood::unordered_flat_map<
        uint64_t,
        robin_hood::unordered_flat_map<
            read_group,
            robin_hood::unordered_flat_map<
            //    std::string, std::pair<seqan::BamAlignmentRecord, uint64_t>>>>
            std::string, std::pair<std::unique_ptr<bam1_t, bam1_t_deleter>, uint64_t>>>>
        read_map;
    robin_hood::unordered_flat_map<
        uint64_t, robin_hood::unordered_flat_map<
                      read_group,
                      robin_hood::unordered_flat_map<std::string, uint64_t>>>
        read_counts;

    auto output_positions = [&read_map, &read_counts, &fun](std::optional<uint64_t> start){
        std::vector<uint64_t> positions;
        positions.reserve(read_map.size());
        for(auto& [k,v] : read_map){
            if(!start.has_value() || k + 1000 < start){
                positions.push_back(k);
            }
        }
        std::sort(positions.begin(), positions.end());
        for(auto p : positions){
            std::vector<read_group> sorted_keys;
            sorted_keys.reserve(read_map.size());
            auto& map_p = read_map[p];
            for(auto& [k,v] : read_map[p]){
                sorted_keys.push_back(k);
            }
            std::sort(sorted_keys.begin(), sorted_keys.end());
            for(auto& k: sorted_keys){
                fun(map_p[k]);
            }
        }
        for(auto p : positions){
            read_map.erase(p);
            read_counts.erase(p);
        }
    };

    //seqan::BamAlignmentRecord record;
    bam1_t* record = bam_init1();
//    while (!seqan::atEnd(file)) {
    while(sam_read1(file, bam_hdr, record) > 0) {
//      seqan::readRecord(record, file);
//      if(!seqan::hasFlagUnmapped(record)) {
      if ((record->core.flag & BAM_FUNMAP) == 0) {
//        cur_ref = record.rID;
        cur_ref = record->core.tid;
        auto* qname = bam_get_qname(record);
        auto umi = get_umi(bam_get_qname(record), ::strlen(qname));
//        auto umi = get_umi(record.qName);
        auto [start, pos, is_spliced] =
            get_read_position(record, opts.soft_clip_threshold);

        if(last_output_pos + 1000 < start || cur_ref != last_ref) {
            output_positions(start);
            last_output_pos = start;
        }

        last_pos = std::max(pos, start);
        last_ref = cur_ref;
        auto key = read_group{
                bam_is_rev(record)
            /*seqan::hasFlagRC(record)*/, opts.spliced && is_spliced != 0,
            static_cast<uint16_t>(opts.read_length ? record->core.l_qseq /*seqan::length(record.seq)*/
                                                   : 0)};
        update_read_map(record, pos, key, std::string{umi}, read_map,
                        read_counts);
      }
    }
    output_positions(std::nullopt);
    bam_destroy1(record);
}


void dedup(const std::string& input, const std::string& output, umi_opts opts) {
  //omp_set_num_threads(1);
  rand_gen.seed(opts.seed);

  samFile* file = hts_open(input.c_str(), "r");
//  seqan::BamFileIn file;
//  if (!seqan::open(file, input.c_str())) {
//    throw std::runtime_error(fmt::format("Could not open file '{}'", input));
//  }

  // read header
  bam_hdr_t* bam_hdr = sam_hdr_read(file);
//  seqan::BamHeader header;
//  seqan::readHeader(header, file);

//  if (seqan::getSortOrder(header) != seqan::BamSortOrder::BAM_SORT_COORDINATE) {
//    throw std::runtime_error(
//        "Deduplication requires a coordinate sorted bam file!");
//  }

  samFile* out = hts_open(output.c_str(), "wb");
  sam_hdr_write(out, bam_hdr);

//  seqan::BamFileOut out(file);
//  if(!seqan::open(out, output.c_str())) {
//      throw std::runtime_error(fmt::format("Could not create file '{}'", output));
//  }
//  seqan::writeHeader(out, header);

  umi_clusterer clusterer(opts.method);
  process_bam_read_chunks(file, bam_hdr, opts, [&clusterer, &out, &bam_hdr](auto& bundle){
      clusterer(bundle, [&out, &bam_hdr](auto& read, auto& umi, auto& count){
        //seqan::writeRecord(out, read);
          sam_write1(out, bam_hdr, read.get());
      });
  });
  bam_hdr_destroy(bam_hdr);

  hts_close(file);
  hts_close(out);
}

}  // namespace fumi_tools
