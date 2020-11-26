#include <fumi_tools/sample_index_map.hpp>
#include <numeric>

#include <fstream>
#include <iostream>

#include <fmt/format.h>
#include <fmt/ostream.h>

#include <ghc/filesystem.hpp>

#include <rapidcsv/rapidcsv.h>

namespace {
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
                                   unsigned int max_errors,
                                   const std::vector<unsigned int>& lanes)
    : max_errors_(max_errors) {
  std::ifstream ifs(sample_sheet);

  if (!ifs) {
    std::cerr << "Failed to open sample sheet '" << sample_sheet
              << "'! Please make sure that it exists and is readable."
              << std::endl;
    std::exit(1);
  }

  // skip until [Data] line, after this we have the sample sheet table
  for (std::string line; std::getline(ifs, line);) {
    if (line.compare(0, 6, "[Data]") == 0) {
      break;
    }
  }

  auto tbl_text = std::string(std::istreambuf_iterator<char>(ifs),
                              std::istreambuf_iterator<char>());
  std::stringstream sstream;
  sstream << tbl_text;
  rapidcsv::Document doc(sstream);
  auto lane_i = doc.GetColumnIdx("Lane");
  if (lane_i == -1 && lanes.empty()) {
    std::cerr << "The lane(s) must be specified either on the command-line or in the sample sheet!"
              << std::endl;
    std::exit(1);
  }
  auto i7_i = doc.GetColumnIdx("index");
  auto i5_i = doc.GetColumnIdx("index2");
  if (i7_i == -1) {
    std::cerr << "Failed to detect necessary i7 index column called 'index' in "
                 "header of provided sample sheet!"
              << std::endl;
    std::exit(1);
  }
  if (i5_i == -1) {
    std::cerr << "Failed to detect necessary i5 index column called 'index'2 "
                 "in header of provided sample sheet!"
              << std::endl;
    std::exit(1);
  }
  auto sample_i = doc.GetColumnIdx("Sample_ID");
  if (sample_i == -1) {
    std::cerr << "Failed to detect necessary sample ID column called "
                 "'Sample_ID' in header of provided sample sheet!"
              << std::endl;
    std::exit(1);
  }

  auto sample_n = doc.GetColumnIdx("Sample_Name");
  if (sample_n == -1) {
    std::cerr << "Failed to detect necessary sample name column called "
                 "'Sample_Name' in header of provided sample sheet!"
              << std::endl;
    std::exit(1);
  }

  auto si = output_pattern.find("%i");
  auto sn = output_pattern.find("%s");
  auto ln = output_pattern.find("%l");
  if (si == nonstd::string_view::npos && sn == nonstd::string_view::npos) {
    std::cerr << "At least the sample name (%s) or sample index (%i) "
                 "placeholders need to be provided in the output"
              << std::endl;
    std::exit(1);
  }

  if (ln == nonstd::string_view::npos) {
    std::cerr << "The lane (%l) placeholder needs to be provided in the output"
              << std::endl;
    std::exit(1);
  }

  auto get_output_filename = [&doc, sample_i, sample_n, output_pattern](
                                 auto i, auto lane) {
    auto output = output_pattern.to_string();
    auto si = output.find("%i");
    if (si != nonstd::string_view::npos) {
      output.replace(si, 2, doc.GetCell<std::string>(sample_i, i));
    }
    auto sn = output.find("%s");
    if (sn != nonstd::string_view::npos) {
      output.replace(sn, 2, doc.GetCell<std::string>(sample_n, i));
    }
    auto ln = output.find("%l");
    if (ln != nonstd::string_view::npos) {
      output.replace(ln, 2, fmt::format("{:03d}", lane));
    }
    return output;
  };

  for (auto i = 0u; i < doc.GetRowCount(); ++i) {
    auto lane = 1u;
    if (lane_i != -1) {
      lane = doc.GetCell<unsigned int>(lane_i, i);
    }
    if (!lanes.empty()) {
      // no lanes in sample sheet but user assigns indices to provided lane(s)
      if (lane_i == -1) {
        for (auto& l : lanes) {
          add_i5_i7_index(doc.GetCell<std::string>(i5_i, i),
                          doc.GetCell<std::string>(i7_i, i), l,
                          get_output_filename(i, l));
        }
        continue;
      } else if (std::find(lanes.begin(), lanes.end(), lane) == lanes.end()) {
        // lanes in sample sheet and user restricts the samples to the provided
        // lane(s)
        continue;
      }
    }
    add_i5_i7_index(doc.GetCell<std::string>(i5_i, i),
                    doc.GetCell<std::string>(i7_i, i), lane,
                    get_output_filename(i, lane));
  }

  // add undetermined files for each lane
  // but only if we have any indices associated with the lane
  for (auto i = 0ul; i < output_files_.size(); ++i) {
    if (!output_files_[i].empty()) {
      auto und_output = output_pattern.to_string();
      auto si = und_output.find("%i");
      if (si != nonstd::string_view::npos) {
        und_output.replace(si, 2, "0");
      }
      auto sn = und_output.find("%s");
      if (sn != nonstd::string_view::npos) {
        und_output.replace(sn, 2, "Undetermined");
      }
      auto ln = und_output.find("%l");
      if (ln != nonstd::string_view::npos) {
        und_output.replace(ln, 2, fmt::format("{:03d}", i + 1));
      }
      output_files_[i].emplace_back(und_output);
    }
  }

  if (i7_indices_.empty()) {
    std::cerr << "Failed to detect any i7 index!" << std::endl;
    std::exit(1);
  }

  i7_length_.resize(i7_indices_.size());
  for (auto i = 0ul; i < i7_indices_.size(); ++i) {
    if (!i7_indices_[i].empty()) {
      i7_length_[i] = i7_indices_[i][0].size();
    }
  }

  for (auto i = 0ul; i < i7_indices_.size(); ++i) {
    for (auto& i7 : i7_indices_[i]) {
      if (i7.size() != i7_length_[i]) {
        std::cerr
            << fmt::format(
                   "Not all i7 indices of lane {:3d} have the sample length!",
                   i + 1)
            << std::endl;
        std::exit(1);
      }
    }
  }

  if (i5_indices_.empty()) {
    std::cerr << "Failed to detect any i5 index!" << std::endl;
    std::exit(1);
  }

  i5_length_.resize(i5_indices_.size());
  for (auto i = 0ul; i < i5_indices_.size(); ++i) {
    if (!i5_indices_[i].empty()) {
      i5_length_[i] = i5_indices_[i][0].size();
    }
  }

  for (auto i = 0ul; i < i5_indices_.size(); ++i) {
    for (auto& i5 : i5_indices_[i]) {
      if (i5.size() != i5_length_[i]) {
        std::cerr
            << fmt::format(
                   "Not all i5 indices of lane {:3d} have the sample length!",
                   i + 1)
            << std::endl;
        std::exit(1);
      }
    }
  }

  // check if we have ambiguous indices when considering mismatches
  for (auto i = 0ul; i < i7_indices_.size(); ++i) {
    for (auto& i7 : i7_indices_[i]) {
      for (auto& other_i7 : i7_indices_[i]) {
        if (i7 != other_i7 &&
            get_num_mismatches(i7, other_i7) <= 2 * max_errors_) {
          std::cerr
              << fmt::format(
                     "Found ambiguous i7 indices in lane {:3d} when allowing "
                     "up to {} mismatches!\ni7 index 1: {}, i7 index 2: {}\n"
                     "Please reduce the number of allowed mismatches.",
                     i + 1, max_errors_, i7, other_i7)
              << std::endl;
          std::exit(1);
        }
      }
    }
  }
  for (auto i = 0ul; i < i5_indices_.size(); ++i) {
    for (auto& i5 : i5_indices_[i]) {
      for (auto& other_i5 : i5_indices_[i]) {
        if (i5 != other_i5 &&
            get_num_mismatches(i5, other_i5) <= 2 * max_errors_) {
          std::cerr
              << fmt::format(
                     "Found ambiguous i5 indices in lane {:3d} when allowing "
                     "up to {} mismatches!\ni5 index 1: {}, i5 index 2: {}\n"
                     "Please reduce the number of allowed mismatches.",
                     i + 1, max_errors_, i5, other_i5)
              << std::endl;
          std::exit(1);
        }
      }
    }
  }
}

uint64_t sample_index_map::find_indices(nonstd::string_view i7,
                                        nonstd::string_view i5,
                                        unsigned int lane) const {
  if (lane > i7_indices_.size() || i7_indices_[lane - 1].empty()) {
    return std::numeric_limits<uint64_t>::max();
  }
  auto it =
      std::find(i7_indices_[lane - 1].begin(), i7_indices_[lane - 1].end(), i7);
  if (it == i7_indices_[lane - 1].end()) {
    auto min_it = i7_indices_[lane - 1].end();
    auto min_mismatches = std::numeric_limits<uint64_t>::max();
    for (auto mit = i7_indices_[lane - 1].begin();
         mit != i7_indices_[lane - 1].end(); ++mit) {
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
  if (it != i7_indices_[lane - 1].end()) {
    auto pos = std::distance(i7_indices_[lane - 1].begin(), it);
    if (i5_indices_[lane - 1][pos] == i5) {
      return pos;
    } else if (get_num_mismatches(i5_indices_[lane - 1][pos], i5) <=
               max_errors_) {
      return pos;
    }
  }
  return output_files_[lane - 1].size() - 1;
}

void sample_index_map::add_i5_i7_index(std::string index5,
                                       std::string index7,
                                       unsigned int lane,
                                       std::string output_filename) {
  if (i7_indices_.size() <= lane - 1) {
    i7_indices_.resize(lane);
  }
  if (i5_indices_.size() <= lane - 1) {
    i5_indices_.resize(lane);
  }
  if (output_files_.size() <= lane - 1) {
    output_files_.resize(lane);
  }
  i7_indices_[lane - 1].push_back(std::move(index7));
  i5_indices_[lane - 1].push_back(std::move(index5));
  output_files_[lane - 1].emplace_back(output_filename);
}

}  // namespace fumi_tools
