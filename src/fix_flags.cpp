#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include <nonstd/string_view.hpp>

#include <cxxopts/cxxopts.hpp>

#include <cpg/cpg.hpp>

#include <fmt/format.h>

#include <gitversion/version.h>

#include <fumi_tools/cast_helper.hpp>
#include <fumi_tools/helper.hpp>

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
      ("input-threads", "", cxxopts::value<uint64_t>()->default_value("1"))
      ("output-threads", "", cxxopts::value<uint64_t>()->default_value("1"))
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
    required_options(
        opts, {"input", "output"});
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
  if(sv == "-"){
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
  namespace  {

  void add_flag(bam1_t* b, uint16_t flag){
      b->core.flag |= flag;
  }

  void remove_flag(bam1_t* b, uint16_t flag){
      b->core.flag &= ~flag;
  }

  void fix_and_output_read_flags(const std::vector<std::unique_ptr<bam1_t, bam1_t_deleter>>& reads, bam_hdr_t* bam_hdr, samFile* outfile){
      // use the mapping quality to determine which alignment should be the primary alignment
      auto best_mapq_i = 0ul;
      auto best_mapq = 0;
      for(auto i = 0ul; i < reads.size(); ++i){
        auto* nh_tag = bam_aux_get(reads[i].get(), "NH");
        auto num_hits = bam_aux2i(nh_tag);
        // adjust number of hits and accordingly number the hit indicies
        if(num_hits != as_signed(reads.size())){
            bam_aux_update_int(reads[i].get(), "NH", as_signed(reads.size()));
            bam_aux_update_int(reads[i].get(), "HI", as_signed(i+1));
        }
        auto qual = reads[i]->core.qual;
        if(qual > best_mapq){
            best_mapq = qual;
            best_mapq_i = i;
        }
      }

      for(auto i = 0ul; i < reads.size(); ++i){
          if(i == best_mapq_i){
              remove_flag(reads[i].get(), BAM_FSECONDARY);
          } else {
              add_flag(reads[i].get(), BAM_FSECONDARY);
          }
          if(sam_write1(outfile, bam_hdr, reads[i].get()) < 0){
              std::cerr << "Failed to write to output file!" << std::endl;
              std::exit(1);
          }
      }
  }

    void process_bam_read_chunks(samFile* infile, bam_hdr_t* bam_hdr, samFile* outfile) {
      std::vector<std::unique_ptr<bam1_t, bam1_t_deleter>> reads;

      cpg::cpg_cfg prog_cfg{};
      prog_cfg.unit = "aln";
      prog_cfg.unit_scale = true;
      prog_cfg.mininterval = 3;
      prog_cfg.desc = "Fix flags";

      auto progress = cpg::cpg(prog_cfg);

      auto last_qname = "";
      bam1_t* record = bam_init1();
      while (sam_read1(infile, bam_hdr, record) > 0) {
        if ((record->core.flag & BAM_FUNMAP) == 0) {
          auto* qname = bam_get_qname(record);
          if(qname != last_qname){
              // fix flags for this chunk of reads and output them
              fix_and_output_read_flags(reads, bam_hdr, outfile);
              reads.clear();
          }
          reads.emplace_back(bam_dup1(record));
        }
        progress.update();
      }
      fix_and_output_read_flags(reads, bam_hdr, outfile);
      bam_destroy1(record);
    }
  }
    void fix_flags(const std::string& input, const std::string& output, uint64_t ithreads, uint64_t othreads){
        samFile* file = hts_open(input.c_str(), "r");

        if (file == nullptr) {
          throw std::runtime_error(fmt::format("Could not open file '{}'", input));
        }

        if(ithreads > 1) {
          hts_set_threads(file, ithreads);
        }

        // read header
        bam_hdr_t* bam_hdr = sam_hdr_read(file);

//        Don't require sorted file, grouping by reads is enough
//        // require sorted bam file
//        if ((bam_hdr->l_text > 3) && (strncmp(bam_hdr->text, "@HD", 3) == 0)) {
//          auto* p = strstr(bam_hdr->text, "\tSO:queryname");
//          auto* q = strchr(bam_hdr->text, '\n');
//          // Looking for SO:queryname within the @HD line only
//          // (e.g. must ignore in a @CO comment line later in header)
//          if ((p == nullptr) || (p >= q)) {
//            throw std::runtime_error(
//                fmt::format("BAM file needs to be read name sorted!"));
//          }
//        } else {
//          throw std::runtime_error(
//              fmt::format("BAM file needs to be read name sorted!"));
//        }

        samFile* out = hts_open(output.c_str(), ends_with(output, ".bam") ? "wb" : "w");
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

        process_bam_read_chunks(file, bam_hdr, out);

        bam_hdr_destroy(bam_hdr);

        hts_close(file);
        hts_close(out);
    }
}


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
                          vm_opts["input-threads"].as<uint64_t>(),
                          vm_opts["output-threads"].as<uint64_t>());
    return 0;
}
