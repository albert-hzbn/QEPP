#pragma once

#include <string>

#include "qe/types.hpp"

namespace qe {

BandData parse_band_table(const std::string& bandPath);

void write_band_plot_bundle(const std::string& bandInputPath,
                            const BandData& bandData,
                            double fermiEv,
                            const std::string& outPrefix);

}  // namespace qe
