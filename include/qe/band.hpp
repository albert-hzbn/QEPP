#pragma once

#include <string>
#include <vector>

#include "qe/types.hpp"

namespace qe {

BandData parse_band_table(const std::string& bandPath);

void write_band_plot_bundle(const std::string& bandInputPath,
                            const BandData& bandData,
                            double fermiEv,
                            const std::string& outPrefix);

// Parse projwfc.x's atomic_proj.xml file.
AtomicProj parse_atomic_proj(const std::string& xmlPath);

// Plot fat band structure: gray band lines + per-group scatter bubbles.
// elements  : filter by element name (e.g. {"Si", "Fe"}). Empty = all elements.
// atomNums  : filter by 1-based atom index. Empty = not used.
// orbitalLs : filter by angular-momentum l (0=s,1=p,2=d,3=f). Empty = all orbitals.
//
// Grouping rules:
//   elements + orbitals  → one group per (element, orbital) pair
//   elements only        → one group per element
//   atoms + orbitals     → one group per (atom, orbital) pair
//   atoms only           → one group per atom
//   orbitals only        → one group per orbital, summed over all atoms
//   none                 → auto: one group per distinct element
void write_fatband_plots(const std::string& bandFilePath,
                         const BandData& bandData,
                         const AtomicProj& proj,
                         double fermiEv,
                         const std::string& outPrefix,
                         const std::vector<std::string>& elements,
                         const std::vector<int>& atomNums,
                         const std::vector<int>& orbitalLs = {});

}  // namespace qe
