#ifndef FUMI_TOOLS_HELPER_HPP
#define FUMI_TOOLS_HELPER_HPP

#include <cstdint>

#include <functional>

#include <htslib/sam.h>

#include <nonstd/string_view.hpp>

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

  read_group() = default;
  read_group(bool is_rev, bool is_spl, int32_t, uint16_t readl)
      :is_reversed(is_rev), is_spliced(is_spl), read_len(readl)
  {}
};

std::ostream& operator<<(std::ostream& out, const read_group& lhs){
  return out << "is_reversed: " << lhs.is_reversed << '\t'
             << "is_splieced: " << lhs.is_spliced << '\t'
             << "read_len: " << lhs.read_len;
}

struct read_group_paired : read_group {
    int32_t template_len;

    read_group_paired() = default;
    read_group_paired(bool is_rev, bool is_spl, int32_t template_l, uint16_t readl)
        :read_group(is_rev, is_spl, template_l, readl), template_len(template_l)
    {}
};

std::ostream& operator<<(std::ostream& out, const read_group_paired& lhs){
  return out << "is_reversed: " << lhs.is_reversed << '\t'
             << "is_splieced: " << lhs.is_spliced << '\t'
             << "read_len: " << lhs.read_len << '\t'
             << "template_len: " << lhs.template_len;
}

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

bool operator<(const read_group_paired& lhs, const read_group_paired& rhs){
    if(lhs.is_reversed != rhs.is_reversed){
        return lhs.is_reversed < rhs.is_reversed;
    } else if(lhs.is_spliced != rhs.is_spliced){
        return lhs.is_spliced < rhs.is_spliced;
    } else if(lhs.template_len != rhs.template_len){
        return lhs.template_len < rhs.template_len;
    } else {
        return lhs.read_len < rhs.read_len;
    }
}

bool operator==(const read_group_paired& lhs, const read_group_paired& rhs) {
  return lhs.is_reversed == rhs.is_reversed &&
         lhs.is_spliced == rhs.is_spliced && lhs.read_len == rhs.read_len &&
         lhs.template_len == rhs.template_len;
}

bool ends_with(const nonstd::string_view str,
               const nonstd::string_view ending) {
  if (str.length() >= ending.length()) {
    return (str.compare(str.length() - ending.length(), ending.length(),
                        ending) == 0);
  } else {
    return false;
  }
}
}
}  // namespace fumi_tools

namespace std {
template <>
struct hash<fumi_tools::read_group> {
  std::size_t operator()(const fumi_tools::read_group& lhs) const noexcept {
    return lhs.is_spliced ^ lhs.is_reversed ^ hash<uint16_t>()(lhs.read_len);
  }
};
template <>
struct hash<fumi_tools::read_group_paired> {
  std::size_t operator()(const fumi_tools::read_group_paired& lhs) const noexcept {
    return lhs.is_spliced ^ lhs.is_reversed ^ hash<uint16_t>()(lhs.read_len) ^ hash<int32_t>()(lhs.template_len);
  }
};

}  // namespace std

#endif  // FUMI_TOOLS_HELPER_HPP
