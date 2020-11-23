#include <fumi_tools/sample_index_map.hpp>
#include <numeric>

#include <fstream>
#include <iostream>

#include <fmt/format.h>
#include <fmt/ostream.h>

#include <ghc/filesystem.hpp>

namespace {
std::vector<std::string> split(const std::string& s, char delimiter) {
  std::vector<std::string> result;
  auto start = 0U;
  auto end = s.find(delimiter);
  while (end != std::string::npos) {
    result.push_back(s.substr(start, end - start));
    start = end + 1;
    end = s.find(',', start);
  }
  result.push_back(s.substr(start, end - start));
  return result;
}

uint64_t get_num_mismatches(nonstd::string_view lhs, nonstd::string_view rhs) {
  auto num_mismatches = 0ul;
  for (auto i = 0ul; i < lhs.size(); ++i) {
    num_mismatches += lhs[i] != rhs[i];
  }
  return num_mismatches;
}
}  // namespace

namespace fumi_tools {
sample_index_map::sample_index_map(const std::string& sample_sheet,
                                   nonstd::string_view output_pattern,
                                   unsigned int max_errors)
    : max_errors_(max_errors) {
  std::ifstream ifs(sample_sheet);

  if (!ifs) {
    std::cerr << "Failed to open sample sheet '" << sample_sheet
              << "'! Please make sure that it exists and is readable."
              << std::endl;
    std::exit(1);
  }

  bool starts = false;
  int64_t sample_i = -1;
  int64_t sample_n = -1;
  int64_t i7_index = -1;
  int64_t i5_index = -1;
  for (std::string line; std::getline(ifs, line);) {
    if (!starts && i7_index == -1 && line.compare(0, 6, "[Data]") == 0) {
      starts = true;
    } else if (starts) {
      auto header = split(line, ',');
      sample_i = std::distance(
          header.begin(), std::find(header.begin(), header.end(), "Sample_ID"));
      sample_n =
          std::distance(header.begin(),
                        std::find(header.begin(), header.end(), "Sample_Name"));
      i7_index = std::distance(
          header.begin(), std::find(header.begin(), header.end(), "index"));
      i5_index = std::distance(
          header.begin(), std::find(header.begin(), header.end(), "index2"));
      starts = false;
    } else if (i7_index != -1) {
      auto els = split(line, ',');
      i7_indices_.push_back(els[i7_index]);
      i5_indices_.push_back(els[i5_index]);
      auto output = output_pattern.to_string();
      auto si = output_pattern.find("%i");
      auto sn = output_pattern.find("%s");
      if (si != nonstd::string_view::npos) {
        output.replace(si, 2, els[sample_i]);
      }
      if (sn != nonstd::string_view::npos) {
        output.replace(sn, 2, els[sample_n]);
      }
      output_files_.push_back(std::make_unique<zstr::ofstream>(output));
      output_filenames_.push_back(output);
    }
  }

  auto und_output = output_pattern.to_string();
  auto si = output_pattern.find("%i");
  auto sn = output_pattern.find("%s");
  if (si != nonstd::string_view::npos) {
    und_output.replace(si, 2, "0");
  }
  if (sn != nonstd::string_view::npos) {
    und_output.replace(sn, 2, "Undetermined");
  }
  output_files_.push_back(std::make_unique<zstr::ofstream>(und_output));
  output_filenames_.push_back(und_output);

  if (i7_indices_.empty()) {
    std::cerr << "Failed to detect any i7 index (no header call 'index' found)!"
              << std::endl;
    std::exit(1);
  }
  i7_length_ = i7_indices_[0].size();

  for (auto& i7 : i7_indices_) {
    if (i7.size() != i7_length_) {
      std::cerr << "Not all i7 indices have the sample length!" << std::endl;
    }
  }

  if (i5_indices_.empty()) {
    std::cerr
        << "Failed to detect any i5 index (no header call 'index2' found)!"
        << std::endl;
    std::exit(1);
  }
  i5_length_ = i5_indices_[0].size();

  for (auto& i5 : i5_indices_) {
    if (i5.size() != i5_length_) {
      std::cerr << "Not all i5 indices have the sample length!" << std::endl;
    }
  }

  // check if we have ambiguous indices when considering mismatches
  for (auto& i7 : i7_indices_) {
    for (auto& other_i7 : i7_indices_) {
      if (i7 != other_i7 &&
          get_num_mismatches(i7, other_i7) <= 2 * max_errors_) {
        std::cerr << fmt::format(
                         "Found ambiguous i7 indices when allowing up to {} "
                         "mismatches!\ni7 index 1: {}, i7 index 2: {}\nPlease "
                         "reduce the number of allowed mismatches.",
                         max_errors_, i7, other_i7)
                  << std::endl;
        std::exit(1);
      }
    }
  }
  for (auto& i5 : i5_indices_) {
    for (auto& other_i5 : i5_indices_) {
      if (i5 != other_i5 &&
          get_num_mismatches(i5, other_i5) <= 2 * max_errors_) {
        std::cerr << fmt::format(
                         "Found ambiguous i5 indices when allowing up to {} "
                         "mismatches!\ni5 index 1: {}, i5 index 2: {}\nPlease "
                         "reduce the number of allowed mismatches.",
                         max_errors_, i5, other_i5)
                  << std::endl;
        std::exit(1);
      }
    }
  }
}

sample_index_map::~sample_index_map() {
  // close all files
  output_files_.clear();
  // delete all empty files
  for (auto& f : output_filenames_) {
    if (ghc::filesystem::is_empty(f)) {
      ghc::filesystem::remove(f);
    }
  }
}

uint64_t sample_index_map::find_indices(nonstd::string_view i7,
                                        nonstd::string_view i5) const {
  auto it = std::find(i7_indices_.begin(), i7_indices_.end(), i7);
  if (it == i7_indices_.end()) {
    auto min_it = i7_indices_.end();
    auto min_mismatches = std::numeric_limits<uint64_t>::max();
    for (auto mit = i7_indices_.begin(); mit != i7_indices_.end(); ++mit) {
      auto nmis = get_num_mismatches(*mit, i7);
      if (nmis < min_mismatches) {
        min_mismatches = nmis;
        min_it = mit;
      }
    }

    if (min_mismatches <= max_errors_) {
      it = min_it;
    }
  }
  if (it != i7_indices_.end()) {
    auto pos = std::distance(i7_indices_.begin(), it);
    if (i5_indices_[pos] == i5) {
      return pos;
    } else if (get_num_mismatches(i5_indices_[pos], i5) <= max_errors_) {
      return pos;
    }
  }
  return output_files_.size() - 1;
}

}  // namespace fumi_tools
