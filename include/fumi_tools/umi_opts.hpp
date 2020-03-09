#ifndef FUMI_TOOLS_UMI_OPTS_HPP
#define FUMI_TOOLS_UMI_OPTS_HPP

#include <string>
#include <vector>

namespace fumi_tools {

struct umi_opts {
  unsigned int max_ham_dist = 1;
  bool read_length = false;
  uint32_t soft_clip_threshold = 4;
  bool spliced = false;
  uint64_t seed = 42;
  std::string method = "unique";
  bool uncompressed = false;
  uint64_t ithreads = 1;
  uint64_t othreads = 1;
  bool paired = false;
  bool ignore_tlen = false;
  std::string unpaired_reads = "use";
  std::string chimeric_pairs = "use";
  std::string unmapped_reads = "discard";
};

}  // namespace fumi_tools

#endif  // FUMI_TOOLS_UMI_OPTS_HPP
