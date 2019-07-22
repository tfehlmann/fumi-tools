#ifndef FUMI_TOOLS_HELPER_HPP
#define FUMI_TOOLS_HELPER_HPP

#include <cstdint>

#include <functional>

#include <htslib/sam.h>

namespace fumi_tools {

namespace {
struct bam1_t_deleter {
    void operator()(bam1_t* lhs) const {
        bam_destroy1(lhs);
    }
};


struct read_group {
  bool is_reversed;
  bool is_spliced;
  uint16_t read_len;
};

bool operator<(const read_group& lhs, const read_group& rhs){
    if(lhs.is_reversed != rhs.is_reversed){
        return lhs.is_reversed < rhs.is_reversed;
    } else if(lhs.is_spliced != rhs.is_spliced){
        return lhs.is_spliced < rhs.is_spliced;
    } else {
        return lhs.read_len < rhs.read_len;
    }
}

bool operator==(const read_group& lhs, const read_group& rhs) {
  return lhs.is_reversed == rhs.is_reversed &&
         lhs.is_spliced == rhs.is_spliced && lhs.read_len == rhs.read_len;
}
}
}  // namespace fumi_tools

namespace std {
template <>
struct hash<fumi_tools::read_group> {
  std::size_t operator()(const fumi_tools::read_group& lhs) const noexcept {
    return lhs.is_spliced ^ lhs.is_spliced ^ hash<uint16_t>()(lhs.read_len);
  }
};

}  // namespace std

#endif  // FUMI_TOOLS_HELPER_HPP
