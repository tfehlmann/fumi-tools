#ifndef FUMI_TOOLS_SAMPLE_INDEX_MAP_HPP
#define FUMI_TOOLS_SAMPLE_INDEX_MAP_HPP

#include <memory>
#include <vector>

#include <nonstd/string_view.hpp>

#include <zstr/zstr.hpp>

namespace fumi_tools {

class zofstream {
 public:
  zofstream(const std::string& filename) : filename_(filename) {}
  zstr::ofstream& operator*() {
    if (strm_ == nullptr) {
      strm_ = std::make_unique<zstr::ofstream>(filename_);
    }
    return *strm_;
  }

  const std::string& get_filename() const;

 private:
  std::unique_ptr<zstr::ofstream> strm_;
  std::string filename_;
};

class sample_index_map {
 public:
  explicit sample_index_map(const std::string& sample_sheet,
                            nonstd::string_view output_pattern,
                            unsigned int max_errors,
                            const std::vector<unsigned int>& lanes);

  uint64_t find_indices(nonstd::string_view i5,
                        nonstd::string_view i7,
                        unsigned int lane = 1) const;
  uint64_t get_i7_length(unsigned int lane = 1) const {
    return i7_length_[lane - 1];
  }

  uint64_t get_i5_length(unsigned int lane = 1) const {
    return i5_length_[lane - 1];
  }

  zstr::ofstream& get_output_file(unsigned int lane, unsigned int pos) const {
    return *output_files_[lane - 1][pos];
  }

 private:
  void add_i5_i7_index(std::string index5,
                       std::string index7,
                       unsigned int lane,
                       std::string output_filename);

  std::vector<std::vector<std::string>> i5_indices_;
  std::vector<std::vector<std::string>> i7_indices_;
  mutable std::vector<std::vector<zofstream>> output_files_;
  std::vector<uint64_t> i7_length_;
  std::vector<uint64_t> i5_length_;
  unsigned int max_errors_;
};

}  // namespace fumi_tools

#endif  // FUMI_TOOLS_SAMPLE_INDEX_MAP_HPP
