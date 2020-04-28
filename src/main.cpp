#include <iostream>

#include <nonstd/string_view.hpp>

#include <cxxopts/cxxopts.hpp>

#include <fmt/format.h>

#include <fumi_tools/version.hpp>

#include <fumi_tools/umi_opts.hpp>
#include <fumi_tools/dedup.hpp>
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

void check_valid_values(const std::string& val, const std::initializer_list<std::string>& valid_values, const std::string& option_name){
  bool valid = std::find(valid_values.begin(), valid_values.end(), val) != valid_values.end();
  if(!valid){
      throw std::runtime_error(fmt::format("Unknown value passed to option '{}': '{}'. Valid options are: {}", option_name, val, fmt::join(valid_values, "|")));
   }
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

auto parse_options(int argc, char* argv[]) {
  cxxopts::Options opts("fumi_tools", "Options");

  fumi_tools::umi_opts umi_opts;

  // clang-format off
  opts.add_options("help")
      ("i,input", "Input SAM or BAM file.", cxxopts::value<std::string>())
      ("o,output", "Output SAM or BAM file. To output SAM on stdout use '-'.", cxxopts::value<std::string>())
//      ("method", "Which method to use to collapse the UMIs. ", cxxopts::value<std::string>()->default_value("")) only unique supported for now
      ("start-only", "Reads only need the same start position and the same UMI to be considered duplicates.")
      ("paired", "Specifiy this option if your alignment file contains paired end reads.")
      ("chimeric-pairs", "How to handle chimeric read pairs. (discard|use)", cxxopts::value<std::string>(umi_opts.chimeric_pairs)->default_value("use"))
      ("unpaired-reads", "How to handle unpaired reads (e.g. mate did not align) (discard|use)", cxxopts::value<std::string>(umi_opts.unpaired_reads)->default_value("use"))
      ("uncompressed", "Output uncompressed BAM.")
      ("seed", "Random number generator seed.", cxxopts::value<uint64_t>(umi_opts.seed)->default_value("42"))
      ("version", "Display version number.")
      ("h,help", "Show this dialog.")
//      ("max-hamming-dist", "Maximum hamming distance for which to collapse umis.", cxxopts::value<uint32_t>(umi_opts.max_ham_dist)->default_value("1")) not yet supported
      ;

  opts.add_options("invisible")
      ("parse_opts", "Only parse options but don't do anything")
      ("input-threads", "", cxxopts::value<uint64_t>(umi_opts.ithreads)->default_value("1"))
      ("output-threads", "", cxxopts::value<uint64_t>(umi_opts.othreads)->default_value("1"))
      ;
  // clang-format on

  try {
    auto copy_argc = argc;
    opts.parse_positional("input");
    opts.parse(copy_argc, argv);
    if (opts["help"].as<bool>()) {
      std::cout << opts.help({"help"}) << std::endl;
      std::exit(0);
    }
    required_options(
        opts, {"input", "output"});
    check_valid_values(umi_opts.chimeric_pairs, {"use", "discard"}, "chimeric-pairs");
    check_valid_values(umi_opts.unpaired_reads, {"use", "discard"}, "unpaired-reads");

    umi_opts.read_length = opts["paired"].as<bool>() ? false : !opts["start-only"].as<bool>();
    umi_opts.uncompressed = opts["uncompressed"].as<bool>();
    umi_opts.paired = opts["paired"].as<bool>();
  } catch (const std::exception& e) {
    if (opts["help"].as<bool>() || argc == 1) {
      std::cout << opts.help({"help"}) << std::endl;
      std::exit(0);
    } else if (opts["version"].as<bool>()) {
      std::cout << "fumi_tools: " << version::VERSION_STRING << std::endl;
      std::exit(0);
    } else {
      std::cout << e.what() << std::endl;
      std::exit(1);
    }
  }

  auto fmt_in = check_format(opts["input"].as<std::string>());
  if (fmt_in == Format::UNKNOWN) {
    std::cerr << "Unknown input format! Needs to be either bam|sam."
              << std::endl;
    std::exit(1);
  }
  auto fmt_out = check_format(opts["output"].as<std::string>());
  if (fmt_out == Format::UNKNOWN) {
    std::cerr << "Unknown output format! Needs to be either bam|sam."
              << std::endl;
    std::exit(1);
  }

  return std::make_pair(opts, umi_opts);
}


}  // namespace

int main(int argc, char* argv[]) {
  // no need to sync
  std::ios_base::sync_with_stdio(false);
  auto vm_opts = parse_options(argc, argv);

  // only parse options and return
  if(vm_opts.first.count("parse_opts") != 0){
      return 0;
  }

  fumi_tools::dedup(
      vm_opts.first["input"].as<std::string>(),
      vm_opts.first["output"].as<std::string>(),
      vm_opts.second);
  return 0;
}
