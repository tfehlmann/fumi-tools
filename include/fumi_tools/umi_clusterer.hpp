#ifndef FUMI_TOOLS_UMI_CLUSTERER_HPP
#define FUMI_TOOLS_UMI_CLUSTERER_HPP

#include <memory>
#include <string>

#include <nonstd/string_view.hpp>

#include <htslib/sam.h>

#include <robin_hood/robin_hood.h>

#include <fumi_tools/helper.hpp>

namespace fumi_tools {
class umi_clusterer {
 public:
  explicit umi_clusterer(nonstd::string_view method = "unique"):
        method_(method) {}

  template <class Fun>
  void operator()(
      const robin_hood::unordered_flat_map<
          std::string, std::pair<std::unique_ptr<bam1_t, bam1_t_deleter>, uint64_t>>&
          bundle,
          Fun fun) {
    if (method_ == "unique") {
        for(auto& umi_info: bundle){
            fun(umi_info.second.first, umi_info.first, umi_info.second.second);
        }
    }
  }

 private:
  std::string method_;
  uint64_t max_umis_per_position_ = 0;
  uint64_t total_umis_per_position_ = 0;
  uint64_t positions_ = 0;
};
}  // namespace fumi_tools

#endif  // FUMI_TOOLS_UMI_CLUSTERER_HPP
