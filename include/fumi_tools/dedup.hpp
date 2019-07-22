#ifndef FUMI_TOOLS_DEDUP_HPP
#define FUMI_TOOLS_DEDUP_HPP

#include <string>

#include <fumi_tools/umi_opts.hpp>

namespace fumi_tools {
    void dedup(const std::string& input, const std::string& output, umi_opts opts);
}


#endif // FUMI_TOOLS_DEDUP_HPP
