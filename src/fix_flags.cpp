#include <fmt/format.h>

#include <cpg/cpg.hpp>
#include <cxxopts/cxxopts.hpp>
#include <fumi_tools/cast_helper.hpp>
#include <fumi_tools/helper.hpp>
#include <fumi_tools/version.hpp>
#include <iostream>
#include <memory>
#include <nonstd/string_view.hpp>
#include <string>
#include <vector>

namespace {

void required_options(cxxopts::Options& opts,
                      std::initializer_list<std::string> req) {
  for (auto& o : req) {
    if (opts.count(o) == 0) {
      throw std::runtime_error(fmt::format("Option '{}' is required!", o));
    }
  }
}

auto parse_options(int argc, char* argv[]) {
  cxxopts::Options opts("fumi_tools", "Options");

  // clang-format off
  opts.add_options()
      ("i,input", "Input SAM or BAM file.", cxxopts::value<std::string>())
      ("o,output", "Output SAM or BAM file.", cxxopts::value<std::string>())
      ("input-threads", "Number of threads to decompress input.", cxxopts::value<uint64_t>()->default_value("1"))
      ("output-threads", "Number of threads to compress output.", cxxopts::value<uint64_t>()->default_value("1"))
      ("sort-adjacent-pairs", "Keep name sorting, but sort pairs such that R2 always follows R1.")
      ("version", "Display version number.")
      ("help", "Show this dialog.")
      ;
  // clang-format on

  try {
    auto copy_argc = argc;
    opts.parse_positional("input");
    opts.parse(copy_argc, argv);
    if (opts["help"].as<bool>()) {
      std::cout << opts.help() << std::endl;
      std::exit(0);
    }
    required_options(opts, {"input", "output"});
  } catch (const std::exception& e) {
    if (opts["help"].as<bool>() || argc == 1) {
      std::cout << opts.help() << std::endl;
      std::exit(0);
    } else if (opts["version"].as<bool>()) {
      std::cout << "fumi_tools: " << version::VERSION_STRING << std::endl;
      std::exit(0);
    } else {
      std::cout << e.what() << std::endl;
      std::exit(1);
    }
  }

  return opts;
}

enum class Format { BAM, SAM, UNKNOWN };

Format check_format(nonstd::string_view sv) {
  if (sv == "-") {
    return Format::SAM;
  } else if (fumi_tools::ends_with(sv, ".bam")) {
    return Format::BAM;
  } else if (fumi_tools::ends_with(sv, ".sam")) {
    return Format::SAM;
  } else {
    return Format::UNKNOWN;
  }
}

}  // namespace

namespace fumi_tools {
namespace {

void add_flag(bam1_t* b, uint16_t flag) {
  b->core.flag |= flag;
}

void remove_flag(bam1_t* b, uint16_t flag) {
  b->core.flag &= ~flag;
}

std::array<int64_t, 4> get_aln_props_wo_flag_info(const bam1_t* lhs,
                                                  bool has_hi) {
  if (lhs == nullptr) {
    return {};
  }
  auto lhs_is_r1 = (lhs->core.flag & BAM_FREAD1) != 0;
  auto lhs_is_r2 = (lhs->core.flag & BAM_FREAD2) != 0;
  auto lhs_hi = 0;
  auto lhs_pos =
      lhs_is_r1 ? lhs->core.pos : lhs_is_r2 ? lhs->core.mpos : lhs->core.pos;
  auto lhs_tid =
      lhs_is_r1 ? lhs->core.tid : lhs_is_r2 ? lhs->core.mtid : lhs->core.tid;
  auto lhs_isize = lhs_is_r1 ? lhs->core.isize
                             : lhs_is_r2 ? -lhs->core.isize : lhs->core.isize;
  if (has_hi) {
    auto* lhs_hi_tag = bam_aux_get(lhs, "HI");
    lhs_hi = bam_aux2i(lhs_hi_tag);
  }
  return {lhs_tid, lhs_pos, lhs_isize, lhs_hi};
}

std::tuple<bool, bool, int64_t, int64_t, int64_t, int64_t> get_aln_props(
    const bam1_t* lhs,
    bool has_hi) {
  if (lhs == nullptr) {
    return {};
  }
  auto lhs_is_r1 = (lhs->core.flag & BAM_FREAD1) != 0;
  auto lhs_is_r2 = (lhs->core.flag & BAM_FREAD2) != 0;
  auto props = get_aln_props_wo_flag_info(lhs, has_hi);
  return std::make_tuple(!lhs_is_r1, !lhs_is_r2, props[0], props[1], props[2],
                         props[3]);
}

struct mapq_stats {
  std::size_t num_r1_reads = 0;
  std::size_t num_r2_reads = 0;
  std::size_t num_other_reads = 0;
  int32_t best_r1_i = -1;
  int32_t best_r2_i = -1;
  int32_t best_other_i = -1;
};

mapq_stats get_best_mapq(
    const std::vector<std::unique_ptr<bam1_t, bam1_t_deleter>>& reads) {
  mapq_stats result{};

  auto best_mapq_r1 = std::numeric_limits<int32_t>::min();
  auto best_mapq_r2 = std::numeric_limits<int32_t>::min();
  auto best_mapq_other = std::numeric_limits<int32_t>::min();

  for (auto i = 0ul; i < reads.size(); ++i) {
    auto qual = reads[i]->core.qual;
    if ((reads[i]->core.flag & BAM_FREAD1) != 0) {
      ++result.num_r1_reads;
      if (qual > best_mapq_r1) {
        best_mapq_r1 = qual;
        result.best_r1_i = i;
      }
    } else if ((reads[i]->core.flag & BAM_FREAD2) != 0) {
      ++result.num_r2_reads;
      if (qual > best_mapq_r2) {
        best_mapq_r2 = qual;
        result.best_r2_i = i;
      }
    } else {
      ++result.num_other_reads;
      if (qual > best_mapq_other) {
        best_mapq_other = qual;
        result.best_other_i = i;
      }
    }
  }
  return result;
}

std::array<int64_t, 3> get_second_best_as(
    const std::vector<std::unique_ptr<bam1_t, bam1_t_deleter>>& reads,
    const mapq_stats& stats) {
  auto has_as = bam_aux_get(reads[0].get(), "AS") != nullptr;
  auto second_best_as_r1 = [&]() {
    if (stats.num_r1_reads == 0 || !has_as) {
      return -1l;
    }
    auto idx = stats.num_r1_reads == 1 ? 0ul : 1ul;
    auto r = reads[idx].get();
    auto r_as_tag = bam_aux_get(r, "AS");
    return bam_aux2i(r_as_tag);
  }();
  auto second_best_as_r2 = [&]() {
    if (stats.num_r2_reads == 0 || !has_as) {
      return -1l;
    }
    auto idx =
        stats.num_r2_reads == 1 ? stats.num_r1_reads : stats.num_r1_reads + 1;
    auto r = reads[idx].get();
    auto r_as_tag = bam_aux_get(r, "AS");
    return bam_aux2i(r_as_tag);
  }();
  auto second_best_as_other = [&]() {
    if (stats.num_other_reads == 0 || !has_as) {
      return -1l;
    }
    auto idx = stats.num_other_reads == 1
                   ? stats.num_r1_reads + stats.num_r2_reads
                   : stats.num_r1_reads + stats.num_r2_reads + 1;
    auto r = reads[idx].get();
    auto r_as_tag = bam_aux_get(r, "AS");
    return bam_aux2i(r_as_tag);
  }();
  return {second_best_as_r1, second_best_as_r2, second_best_as_other};
}

std::size_t set_primary_alignment(
    const std::vector<std::unique_ptr<bam1_t, bam1_t_deleter>>& reads,
    const mapq_stats& stats,
    bool has_hi) {
  std::size_t distinct_alignments = 0;
  std::size_t r1_idx = 0;
  auto r2_idx = stats.num_r1_reads;
  while (r1_idx < stats.num_r1_reads ||
         r2_idx < stats.num_r1_reads + stats.num_r2_reads) {
    bam1_t* r1 = nullptr;
    bam1_t* r2 = nullptr;
    if (r1_idx < stats.num_r1_reads) {
      r1 = reads[r1_idx].get();
    }

    if (r2_idx < stats.num_r1_reads + stats.num_r2_reads) {
      r2 = reads[r2_idx].get();
    }
    auto r1_props = get_aln_props_wo_flag_info(r1, has_hi);
    auto r2_props = get_aln_props_wo_flag_info(r2, has_hi);

    if (r1 != nullptr && (r2 == nullptr || r1_props < r2_props)) {
      ++r1_idx;
      ++distinct_alignments;
      if (stats.best_r1_i >= 0) {
        if (r1_idx == as_unsigned(stats.best_r1_i)) {
          remove_flag(r1, BAM_FSECONDARY);
        } else {
          add_flag(r1, BAM_FSECONDARY);
        }
      }
    } else if (r2 != nullptr && (r1 == nullptr || r1_props > r2_props)) {
      ++distinct_alignments;
      ++r2_idx;
      add_flag(r2, BAM_FSECONDARY);
    } else {
      ++distinct_alignments;
      if (r1_idx == as_unsigned(stats.best_r1_i)) {
        remove_flag(r1, BAM_FSECONDARY);
        remove_flag(r2, BAM_FSECONDARY);
      } else {
        add_flag(r1, BAM_FSECONDARY);
        add_flag(r2, BAM_FSECONDARY);
      }
      ++r1_idx;
      ++r2_idx;
    }
  }

  // if we don't have any r1, r2 reads there can only be 1 primary alignment
  for (std::size_t i = stats.num_r1_reads + stats.num_r2_reads;
       i < reads.size(); ++i) {
    if (stats.best_r1_i == -1 && stats.best_r2_i == -1 &&
        stats.best_other_i >= 0) {
      if (i == as_unsigned(stats.best_other_i)) {
        remove_flag(reads[i].get(), BAM_FSECONDARY);
      } else {
        add_flag(reads[i].get(), BAM_FSECONDARY);
      }
    }
  }

  distinct_alignments += stats.num_other_reads;
  return distinct_alignments;
}

void update_xs_nh_hi_fields(
    const std::vector<std::unique_ptr<bam1_t, bam1_t_deleter>>& reads,
    const mapq_stats& stats,
    std::size_t distinct_alignments,
    std::array<int64_t, 3>& second_best_as) {
  auto has_xs = bam_aux_get(reads[0].get(), "XS") != nullptr;
  auto has_as = bam_aux_get(reads[0].get(), "AS") != nullptr;
  auto has_nh = bam_aux_get(reads[0].get(), "NH") != nullptr;
  // update second best alignment score (XS), number of hits (NH) and hit index
  // (HI) fields
  for (auto i = 0ul; i < reads.size(); ++i) {
    if (has_nh) {
      if ((reads[i]->core.flag & BAM_FREAD1) != 0) {
        bam_aux_update_int(reads[i].get(), "NH",
                           as_signed(distinct_alignments));
        bam_aux_update_int(reads[i].get(), "HI", as_signed(i + 1));
      } else if ((reads[i]->core.flag & BAM_FREAD2) != 0) {
        bam_aux_update_int(reads[i].get(), "NH",
                           as_signed(distinct_alignments));
        bam_aux_update_int(reads[i].get(), "HI",
                           as_signed(i + 1 - stats.num_r1_reads));
      } else {
        bam_aux_update_int(reads[i].get(), "NH",
                           as_signed(distinct_alignments));
        bam_aux_update_int(
            reads[i].get(), "HI",
            as_signed(i + 1 - stats.num_r1_reads - stats.num_r2_reads));
      }
    }
    if (has_as && has_xs) {
      if ((reads[i]->core.flag & BAM_FREAD1) != 0) {
        bam_aux_update_int(reads[i].get(), "XS", second_best_as[0]);
      } else if ((reads[i]->core.flag & BAM_FREAD2) != 0) {
        bam_aux_update_int(reads[i].get(), "XS", second_best_as[1]);
      } else {
        bam_aux_update_int(reads[i].get(), "XS", second_best_as[2]);
      }
    }
  }
}

bool get_pattern_code(uint32_t flag) {
  if ((flag & BAM_FREAD1) != 0)
    return (flag & BAM_FREVERSE) != 0;
  else
    return (flag & BAM_FREVERSE) == 0;
}

void fix_and_output_read_flags(
    std::vector<std::unique_ptr<bam1_t, bam1_t_deleter>>& reads,
    bam_hdr_t* bam_hdr,
    bool rsem_sort,
    samFile* outfile) {
  if (reads.empty()) {
    return;
  }

  // use the mapping quality to determine which alignment should be the primary
  // alignment
  // collect second best alignment score, if available, to update the XS tag
  // accordingly
  auto has_hi = bam_aux_get(reads[0].get(), "HI") != nullptr;

  // quick path
  if (reads.size() == 1) {
    auto r = reads[0].get();
    // set primary aln
    remove_flag(r, BAM_FSECONDARY);
    auto aux_xs = bam_aux_get(r, "XS");
    auto has_xs = aux_xs != nullptr;
    auto aux_as = bam_aux_get(r, "AS");
    auto has_as = aux_as != nullptr;
    auto aux_nh = bam_aux_get(r, "NH");
    auto has_nh = aux_nh != nullptr;

    if (has_as && has_xs) {
      // second best aln score is self
      auto update = bam_aux2i(aux_as);
      if (update != bam_aux2i(aux_xs)) {
        bam_aux_update_int(r, "XS", bam_aux2i(aux_as));
      }
    }
    if (has_nh) {
      if (bam_aux2i(aux_nh) != 1) {
        bam_aux_update_int(r, "NH", 1);
      }
      auto aux_hi = bam_aux_get(r, "HI");

      if (aux_hi != nullptr && bam_aux2i(aux_hi) != 1) {
        bam_aux_update_int(r, "HI", 1);
      }
    }
    if (sam_write1(outfile, bam_hdr, r) < 0) {
      std::cerr << "Failed to write to output file!" << std::endl;
      std::exit(1);
    }
    return;
  }

  // order so that we have first r1, then r2, then unpaired
  // r1 and r2 ordered the same way such that if they are paired they come in
  // the same order possibility 1: only r1 without paired mate 2: only r2
  // without paired mate 3: r1 & r2 paired 4: not paired
  std::stable_sort(reads.begin(), reads.end(), [has_hi](auto& lhs, auto& rhs) {
    return get_aln_props(lhs.get(), has_hi) < get_aln_props(rhs.get(), has_hi);
  });

  // get best mapping qualities for r1, r2 and rest
  auto mapq_stats = get_best_mapq(reads);
  auto second_best_as = get_second_best_as(reads, mapq_stats);

  auto distinct_alignments = set_primary_alignment(reads, mapq_stats, has_hi);
  update_xs_nh_hi_fields(reads, mapq_stats, distinct_alignments,
                         second_best_as);

  if (rsem_sort) {
    auto rsem_less = [](const auto& lhs, const auto& rhs) {
      auto lhsp = std::minmax(lhs->core.pos, lhs->core.mpos);
      auto rhsp = std::minmax(rhs->core.pos, rhs->core.mpos);
      auto lhspat = get_pattern_code(lhs->core.flag);
      auto rhspat = get_pattern_code(rhs->core.flag);

      if (lhs->core.tid != rhs->core.tid) {
        return lhs->core.tid < rhs->core.tid;
      }
      if (lhsp.first != rhsp.first) {
        return lhsp.first < rhsp.first;
      }
      if (lhsp.second != rhsp.second) {
        return lhsp.second < rhsp.second;
      }
      return lhspat < rhspat;
    };
    std::sort(reads.begin(), reads.end(), rsem_less);
  } else {
    auto samtools_less = [](const auto& lhs, const auto& rhs) {
      return (lhs->core.flag & 0xc0) < (rhs->core.flag & 0xc0);
    };
    std::sort(reads.begin(), reads.end(), samtools_less);
  }
  for (auto& r : reads) {
    if (sam_write1(outfile, bam_hdr, r.get()) < 0) {
      std::cerr << "Failed to write to output file!" << std::endl;
      std::exit(1);
    }
  }
}

nonstd::string_view get_canonical_name(const bam1_t* record) {
  // keep only up to first whitespace in case the aligner did not trim the end
  // off
  nonstd::string_view qname(bam_get_qname(record),
                            std::strlen(bam_get_qname(record)));
  auto end = qname.find_last_not_of(" \t\n");
  if (end != nonstd::string_view::npos) {
    qname.remove_suffix(qname.size() - end - 1);
  }
  return qname;
}

void process_bam_read_chunks(samFile* infile,
                             bam_hdr_t* bam_hdr,
                             bool sort_rsem,
                             samFile* outfile) {
  std::vector<std::unique_ptr<bam1_t, bam1_t_deleter>> reads;

  cpg::cpg_cfg prog_cfg{};
  prog_cfg.unit = "aln";
  prog_cfg.unit_scale = true;
  prog_cfg.mininterval = 3;
  prog_cfg.desc = "Fix flags";

  auto progress = cpg::cpg(prog_cfg);

  std::string last_qname;
  bam1_t* record = bam_init1();
  while (sam_read1(infile, bam_hdr, record) > 0) {
    if ((record->core.flag & BAM_FUNMAP) == 0) {
      auto qname = get_canonical_name(record);
      if (qname != last_qname) {
        // fix flags for this chunk of reads and output them
        fix_and_output_read_flags(reads, bam_hdr, sort_rsem, outfile);
        reads.clear();
        last_qname.assign(qname.begin(), qname.end());
      }
      reads.emplace_back(bam_dup1(record));
    }
    progress.update();
  }
  fix_and_output_read_flags(reads, bam_hdr, sort_rsem, outfile);
  bam_destroy1(record);
}

void fix_flags(const std::string& input,
               const std::string& output,
               bool sort_rsem,
               uint64_t ithreads,
               uint64_t othreads) {
  samFile* file = hts_open(input.c_str(), "r");

  if (file == nullptr) {
    throw std::runtime_error(fmt::format("Could not open file '{}'", input));
  }

  if (ithreads > 1) {
    hts_set_threads(file, ithreads);
  }

  // read header
  bam_hdr_t* bam_hdr = sam_hdr_read(file);

  samFile* out =
      hts_open(output.c_str(), ends_with(output, ".bam") ? "wb" : "w");
  if (out == nullptr) {
    throw std::runtime_error(fmt::format("Could not open file '{}'", output));
  }

  if (othreads > 1) {
    hts_set_threads(out, othreads);
  }

  if (sam_hdr_write(out, bam_hdr) < -1) {
    throw std::runtime_error(
        fmt::format("Could not write header to file '{}'", output));
  }

  process_bam_read_chunks(file, bam_hdr, sort_rsem, out);

  bam_hdr_destroy(bam_hdr);

  hts_close(file);
  hts_close(out);
}
}  // namespace
}  // namespace fumi_tools

int main(int argc, char* argv[]) {
  // no need to sync
  std::ios_base::sync_with_stdio(false);
  auto vm_opts = parse_options(argc, argv);
  auto fmt_in = check_format(vm_opts["input"].as<std::string>());
  if (fmt_in == Format::UNKNOWN) {
    std::cerr << "Unknown input format! Needs to be either bam|sam."
              << std::endl;
    return 1;
  }
  auto fmt_out = check_format(vm_opts["output"].as<std::string>());
  if (fmt_out == Format::UNKNOWN) {
    std::cerr << "Unknown output format! Needs to be either bam|sam."
              << std::endl;
    return 1;
  }
  fumi_tools::fix_flags(vm_opts["input"].as<std::string>(),
                        vm_opts["output"].as<std::string>(),
                        vm_opts["sort-adjacent-pairs"].as<bool>(),
                        vm_opts["input-threads"].as<uint64_t>(),
                        vm_opts["output-threads"].as<uint64_t>());
  return 0;
}
