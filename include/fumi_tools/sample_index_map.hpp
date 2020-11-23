#ifndef FUMI_TOOLS_SAMPLE_INDEX_MAP_HPP
#define FUMI_TOOLS_SAMPLE_INDEX_MAP_HPP

#include <vector>

#include <nonstd/string_view.hpp>

#include <zstr/zstr.hpp>

namespace fumi_tools {
class sample_index_map {
 public:
  explicit sample_index_map(const std::string& sample_sheet,
                            nonstd::string_view output_pattern,
                            unsigned int max_errors);

  ~sample_index_map();

  uint64_t find_indices(nonstd::string_view i5, nonstd::string_view i7) const;
  uint64_t get_i7_length() const { return i7_length_; }
  uint64_t get_i5_length() const { return i5_length_; }
  const std::vector<std::unique_ptr<zstr::ofstream>>& get_output_files() const {
    return output_files_;
  }

 private:
  std::vector<std::string> i5_indices_;
  std::vector<std::string> i7_indices_;
  std::vector<std::unique_ptr<zstr::ofstream>> output_files_;
  std::vector<std::string> output_filenames_;
  uint64_t i7_length_;
  uint64_t i5_length_;
  unsigned int max_errors_;
};

}  // namespace fumi_tools

#endif  // FUMI_TOOLS_SAMPLE_INDEX_MAP_HPP
