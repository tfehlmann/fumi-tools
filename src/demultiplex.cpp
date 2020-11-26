#include <condition_variable>
#include <cpg/cpg.hpp>
#include <cxxopts/cxxopts.hpp>
#include <fumi_tools/cast_helper.hpp>
#include <fumi_tools/version.hpp>
#include <iostream>
#include <memory>
#include <mutex>
#include <nonstd/string_view.hpp>
#include <queue>
#include <string>
#include <thread>
#include <vector>

//#include <tbb/tbb.h>

#include <zstr/zstr.hpp>

#include <fmt/format.h>
#include <fmt/ostream.h>

#include <fumi_tools/sample_index_map.hpp>

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
      ("i,input", "Input FASTQ file.", cxxopts::value<std::string>())
      ("s,sample-sheet", "Sample Sheet in Illumina format", cxxopts::value<std::string>())
      ("o,output", "Output FASTQ file pattern, optionally gzip compressed. Use %i as placeholder for the sample index specified in the sample sheet and %s for the sample name (e.g. demultiplexed_reads/%i_S%s_L001_R1.fastq.gz).", cxxopts::value<std::string>())
      ("e,max-errors", "Maximum allowed number of errors (mismatches per default).", cxxopts::value<unsigned int>()->default_value("1"))
      ("format-umi", "Add UMI to the end of the FASTQ header, as expected by fumi_tools dedup")
      ("l,lane", "Optionally specify on which lane the samples provided in the sample sheet ran. Can be specified multiple times to pass several lanes. This option takes precedence on the Lane column of the sample sheet.", cxxopts::value<std::vector<unsigned int>>())
      ("threads", "Number of threads.", cxxopts::value<unsigned int>()->default_value("1"))
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
    required_options(opts, {"input", "sample-sheet", "output"});
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

bool check_format(nonstd::string_view sv) {
  return sv.ends_with(".fastq.gz") || sv.ends_with(".fq.gz") ||
         sv.ends_with(".fastq") || sv.ends_with(".fq");
}

}  // namespace

namespace fumi_tools {
namespace {

struct fq_entry {
  std::string header;
  std::string seq;
  std::string desc;
  std::string qual;
};

std::ostream& operator<<(std::ostream& os, const fq_entry& e) {
  return os << e.header << '\n'
            << e.seq << '\n'
            << e.desc << '\n'
            << e.qual << '\n';
}

// void demultiplex_serial(const std::string& input, const sample_index_map&
// map) {
//  zstr::ifstream ifs(input);

//  const auto& output_files = map.get_output_files();
//  uint64_t i = 0;
//  fq_entry current{};
//  for (std::string line; std::getline(ifs, line); ++i) {
//    if (i % 4 == 0) {
//      current.header = line;
//    } else if (i % 4 == 1) {
//      current.seq = line;
//    } else if (i % 4 == 2) {
//      current.desc = line;
//    } else if (i % 4 == 3) {
//      current.qual = line;
//      auto i7_start = current.header.rfind(":");
//      if (i7_start != std::string::npos) {
//        i7_start++;
//      }
//      auto i7 = nonstd::string_view(current.header.c_str() + i7_start,
//                                    map.get_i7_length());
//      auto i5 = nonstd::string_view(
//          current.header.c_str() + current.header.size() -
//          map.get_i5_length(), map.get_i5_length());
//      auto pos = map.find_indices(i7, i5);
//      if (pos != std::numeric_limits<uint64_t>::max()) {
//        (*output_files[pos]) << current;
//      }
//    }
//  }
//}

// void demultiplex_parallel(const std::string& input,
//                          const sample_index_map& map) {
//  zstr::ifstream ifs(input);

//  const auto& output_files = map.get_output_files();
//  std::vector<tbb::spin_mutex> mutexes(output_files.size());

//  std::vector<fq_entry> fq_entries;
//  const uint64_t fq_entries_chunk_size = 1024;
//  fq_entries.reserve(fq_entries_chunk_size);
//  tbb::parallel_pipeline(
//      128,
//      tbb::make_filter<void, std::vector<fq_entry>>(
//          tbb::filter::serial,
//          [&fq_entries, &ifs](tbb::flow_control& fc) {
//            std::vector<fq_entry> fq_entries;
//            uint64_t i = 0;
//            fq_entry current{};
//            for (std::string line; std::getline(ifs, line); ++i) {
//              if (i % 4 == 0) {
//                current.header = line;
//              } else if (i % 4 == 1) {
//                current.seq = line;
//              } else if (i % 4 == 2) {
//                current.desc = line;
//              } else if (i % 4 == 3) {
//                current.qual = line;
//                fq_entries.push_back(std::move(current));
//                if (fq_entries.size() == fq_entries_chunk_size) {
//                  return fq_entries;
//                }
//              }
//            }
//            fc.stop();
//            return fq_entries;
//          }) &
//          tbb::make_filter<std::vector<fq_entry>, void>(
//              tbb::filter::parallel,
//              [&map, &output_files, &mutexes](std::vector<fq_entry> entries) {
//                std::vector<std::vector<fq_entry>> out_entries(
//                    output_files.size());
//                for (auto& e : entries) {
//                  auto i7_start = e.header.rfind(":");
//                  if (i7_start != std::string::npos) {
//                    i7_start++;
//                  }
//                  auto i7 = nonstd::string_view(e.header.c_str() + i7_start,
//                                                map.get_i7_length());
//                  auto i5 = nonstd::string_view(
//                      e.header.c_str() + e.header.size() -
//                      map.get_i5_length(), map.get_i5_length());
//                  auto pos = map.find_indices(i7, i5);
//                  if (pos != std::numeric_limits<uint64_t>::max()) {
//                    // out_entries[] e;
//                  }
//                }
//              }));
//}

unsigned int extract_lane(nonstd::string_view header) {
  auto pos = header.find(":");
  pos = header.find(":", pos + 1);
  pos = header.find(":", pos + 1);
  auto lane = std::strtoul(header.data() + pos + 1, nullptr, 10);
  if (lane == 0) {
    throw std::runtime_error(
        fmt::format("Lane could not be extracted from header: {}", header));
  }
  return lane;
}

void demultiplex_parallel2(const std::string& input,
                           const sample_index_map& map,
                           bool format_umi,
                           unsigned int threads) {
  zstr::ifstream ifs(input);

  std::vector<std::queue<std::vector<std::tuple<uint32_t, uint32_t, fq_entry>>>>
      t_queues(threads);
  std::vector<std::mutex> t_mutexes(threads);
  std::vector<std::condition_variable> t_cvs(threads);

  bool is_done = false;
  std::vector<std::thread> out_threads;
  out_threads.reserve(threads);
  for (auto i = 0ul; i < threads; ++i) {
    out_threads.emplace_back(
        [i, &map, &t_mutexes, &t_cvs, &t_queues, &is_done]() {
          while (true) {
            std::unique_lock<std::mutex> _(t_mutexes[i]);
            t_cvs[i].wait(_, [i, &t_queues] { return !t_queues[i].empty(); });
            auto fq_entries = std::move(t_queues[i].front());
            t_queues[i].pop();
            _.unlock();
            for (auto& fq_entry : fq_entries) {
              map.get_output_file(std::get<0>(fq_entry), std::get<1>(fq_entry))
                  << std::get<2>(fq_entry);
            }
            if (is_done && t_queues[i].empty()) {
              break;
            }
          }
        });
  }

  uint64_t i = 0;
  fq_entry current{};
  for (std::string line; std::getline(ifs, line); ++i) {
    if (i % 4 == 0) {
      current.header = line;
    } else if (i % 4 == 1) {
      current.seq = line;
    } else if (i % 4 == 2) {
      current.desc = line;
    } else if (i % 4 == 3) {
      current.qual = line;
      auto i7_start = current.header.rfind(":");
      if (i7_start != std::string::npos) {
        i7_start++;
      } else {
        std::cerr << "Could not find i7 index (no ':' found in line" << line
                  << ")!" << std::endl;
        std::exit(1);
      }
      auto i7 = nonstd::string_view(current.header.c_str() + i7_start,
                                    map.get_i7_length());
      auto i5 = nonstd::string_view(
          current.header.c_str() + current.header.size() - map.get_i5_length(),
          map.get_i5_length());
      auto lane = extract_lane(current.header);
      if (format_umi) {
        auto umi_length = current.header.size() - map.get_i5_length() -
                          i7_start - map.get_i7_length() - 1;
        auto umi = nonstd::string_view(
            current.header.c_str() + i7_start + map.get_i7_length(),
            umi_length);
        current.header += fmt::format("_{}", umi);
      }
      auto pos = map.find_indices(i7, i5, lane);
      if (pos != std::numeric_limits<uint64_t>::max()) {
        std::lock_guard<std::mutex> _(t_mutexes[pos % threads]);
        if (t_queues[pos % threads].empty()) {
          t_queues[pos % threads].push({});
        }
        t_queues[pos % threads].front().push_back(
            std::make_tuple(lane, pos, current));
        if (t_queues[pos % threads].front().size() > 4096) {
          t_cvs[pos % threads].notify_one();
        }
      }
    }
  }
  is_done = true;

  for (auto& cv : t_cvs) {
    cv.notify_one();
  }
  for (auto& t : out_threads) {
    t.join();
  }
}

}  // namespace
}  // namespace fumi_tools

int main(int argc, char* argv[]) {
  // no need to sync
  std::ios_base::sync_with_stdio(false);
  auto vm_opts = parse_options(argc, argv);
  auto in_ok = check_format(vm_opts["input"].as<std::string>());
  if (!in_ok) {
    std::cerr << "Unknown input format! Needs to be either fastq[.gz]|fq[.gz]."
              << std::endl;
    return 1;
  }

  // tbb::task_scheduler_init init(vm_opts["threads"].as<unsigned int>());

  fumi_tools::sample_index_map map(
      vm_opts["sample-sheet"].as<std::string>(),
      vm_opts["output"].as<std::string>(),
      vm_opts["max-errors"].as<unsigned int>(),
      vm_opts["lane"].as<std::vector<unsigned int>>());
  fumi_tools::demultiplex_parallel2(vm_opts["input"].as<std::string>(), map,
                                    vm_opts["format-umi"].as<bool>(),
                                    vm_opts["threads"].as<unsigned int>());
  return 0;
}
