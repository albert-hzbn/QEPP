#pragma once

#include <string>

namespace qe {

// Pre-processing: write a pp.x input for STM simulation (ILDOS, plot_num=5).
//   scfInputPath  — QE SCF input; used to extract prefix, outdir
//   outputDir     — directory where the pp.x input is written (created if absent)
//   biasEv        — bias voltage in eV; positive = empty states [0, bias],
//                   negative = filled states [bias, 0]
void write_stm_pp_input(const std::string& scfInputPath,
                        const std::string& outputDir,
                        double biasEv);

// Post-processing: read the cube produced by pp.x STM (plot_num=5) and plot a
// 2D constant-height XY map.
//   cubePath   — path to the .cube file
//   outPrefix  — output filename prefix
//   heightAng  — tip height in Å above the cube origin; -1 = use midplane z-index
void write_stm_plot(const std::string& cubePath,
                    const std::string& outPrefix,
                    double heightAng = -1.0);

}  // namespace qe
