#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <iostream>
#include <memory>
#include <mutex>
#include <nonstd/string_view.hpp>
#include <queue>
#include <deque>
#include <string>
#include <thread>
#include <vector>

//#include <tbb/tbb.h>

#include <cpg/cpg.hpp>
#include <cxxopts/cxxopts.hpp>

#include <fumi_tools/cast_helper.hpp>
#include <fumi_tools/version.hpp>

#include <zstr/zstr.hpp>

#include <fmt/format.h>
#include <fmt/ostream.h>

#include <fumi_tools/sample_index_map.hpp>

using namespace std::chrono_literals;
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
      ("tag-umi", "Add UMI to the read ID by adding :FUMI|<UMI_SEQ>| instead of a simple underscore.")
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
                           bool tag_umi,
                           unsigned int threads) {
  zstr::ifstream ifs(input);

  std::vector<std::deque<std::vector<std::tuple<uint32_t, uint32_t, std::string>>>>
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
                    t_cvs[i].wait(_, [i, &t_queues, &is_done] {
                        return !t_queues[i].empty() || is_done;
                    });
                    if (is_done && t_queues[i].empty()) {
                      break;
                    }
                    auto fq_entries = std::move(t_queues[i].front());
                    t_queues[i].pop_front();
                    _.unlock();
                    for (auto& fq_entry : fq_entries) {
                      map.get_output_file(std::get<0>(fq_entry),
                      std::get<1>(fq_entry))
                          << std::get<2>(fq_entry);
                    }
                    if (is_done && t_queues[i].empty()) {
                      break;
                    }
                  }
                }
    );
  }

  uint64_t mem_limit_per_thread = 1000 * 1024 * 1024;  // 1 GB
  uint64_t i = 0;
  fq_entry current{};
  std::string entry;
  cpg::cpg_cfg pcfg;
  pcfg.desc = "Demultiplexing";
  pcfg.unit = "reads";
  pcfg.unit_scale = true;
  pcfg.dynamic_ncols = true;
  auto progress = cpg::cpg(pcfg);
  constexpr uint64_t progress_interval = 10000000;
  uint64_t read_count = 0;
  std::vector<uint64_t> skipped_lanes;
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
      auto lane = extract_lane(current.header);
      if (!map.has_lane(lane)) {
        if (lane >= skipped_lanes.size()) {
          skipped_lanes.resize(lane + 1, 0);
        }
        if (skipped_lanes[lane] == 0) {
          fmt::print(stderr,
              "Warning: encountered reads from lane {} which is not in the "
              "sample sheet. These reads will be skipped.\n", lane);
        }
        ++skipped_lanes[lane];
        ++read_count;
        if (read_count % progress_interval == 0) {
          progress.update(progress_interval);
        }
        continue;
      }
      auto i7 = nonstd::string_view(current.header.c_str() + i7_start,
                                    map.get_i7_length(lane));
      auto i5 =
          nonstd::string_view(current.header.c_str() + current.header.size() -
                                  map.get_i5_length(lane),
                              map.get_i5_length(lane));

      if (format_umi) {
        auto umi_length = current.header.size() - map.get_i5_length(lane) -
                          i7_start - map.get_i7_length(lane) - 1;
        auto umi = nonstd::string_view(
            current.header.c_str() + i7_start + map.get_i7_length(lane),
            umi_length);
        if (tag_umi) {
          current.header += fmt::format(":FUMI|{}|", umi);
        } else {
          current.header += fmt::format("_{}", umi);
        }
      }
      auto pos = map.find_indices(i7, i5, lane);
      if (pos != std::numeric_limits<uint64_t>::max()) {
        std::unique_lock<std::mutex> _(t_mutexes[pos % threads]);
        if (t_queues[pos % threads].empty()) {
          t_queues[pos % threads].push_back({});
        }
        auto q_entry = fmt::format("{}\n{}\n{}\n{}\n",
                                 current.header, current.seq, current.desc, current.qual);
        auto q_size = q_entry.size() + sizeof(q_entry);
        t_queues[pos % threads].front().push_back(
            std::make_tuple(lane, pos, std::move(q_entry)));
        if (t_queues[pos % threads].front().size() > 4096) {
          t_cvs[pos % threads].notify_one();
          while(!t_queues[pos % threads].empty() && t_queues[pos % threads].front().size()  * q_size * 2 > mem_limit_per_thread) {
              _.unlock();
              std::this_thread::sleep_for(300ms);
              _.lock();
          }
        }
      }
      ++read_count;
      if (read_count % progress_interval == 0) {
        progress.update(progress_interval);
      }
    }
  }
  progress.update(read_count % progress_interval);
  for (unsigned int lane = 0; lane < skipped_lanes.size(); ++lane) {
    if (skipped_lanes[lane] > 0) {
      fmt::print(stderr, "\nWarning: skipped {} reads from lane {} "
          "(not in the sample sheet)\n", skipped_lanes[lane], lane);
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
  fumi_tools::demultiplex_parallel2(
      vm_opts["input"].as<std::string>(), map, vm_opts["format-umi"].as<bool>(),
      vm_opts["tag-umi"].as<bool>(), vm_opts["threads"].as<unsigned int>());
  return 0;
}
