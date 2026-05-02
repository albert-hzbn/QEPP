#pragma once

#include <string>
#include <vector>

#include "qe/types.hpp"

namespace qe {

// ── I/O ──────────────────────────────────────────────────────────────────────

// Parse a QE bands.x output file (.dat.gnu or 2-column k/E format).
BandData parse_band_table(const std::string& bandPath);

// Parse projwfc.x's atomic_proj.xml file.
AtomicProj parse_atomic_proj(const std::string& xmlPath);

// ── Algorithm ────────────────────────────────────────────────────────────────

// A named group of wavefunction indices used for fat band plotting.
struct FatBandGroup {
    std::string      name;
    std::vector<int> wfcIdx;
};

// Build projection groups from an AtomicProj using optional filters.
//   elements  : group by element name (e.g. {"Si", "Fe"}). Empty = not used.
//   atomNums  : group by 1-based atom index. Empty = not used.
//   orbitalLs : group by angular-momentum l (0=s,1=p,2=d,3=f). Empty = not used.
//
// Grouping rules (filters are combinable):
//   elements + orbitals  → one group per (element, orbital) pair
//   elements only        → one group per element
//   atoms + orbitals     → one group per (atom, orbital) pair
//   atoms only           → one group per atom
//   orbitals only        → one group per orbital, summed over all atoms
//   none                 → auto: one group per distinct element
std::vector<FatBandGroup> build_fatband_groups(
    const AtomicProj& proj,
    const std::vector<std::string>& elements,
    const std::vector<int>& atomNums,
    const std::vector<int>& orbitalLs);

// ── Plotting ─────────────────────────────────────────────────────────────────

// Write a band-structure data file (.dat) and a PNG image.
void write_band_plot_bundle(const std::string& bandInputPath,
                            const BandData& bandData,
                            double fermiEv,
                            const std::string& outPrefix);

// Plot fat (projected) band structure. Gray lines show all bands; each group
// is overlaid in a distinct color. Produces one combined PNG.
void write_fatband_plots(const BandData& bandData,
                         const std::vector<FatBandGroup>& groups,
                         double fermiEv,
                         const std::string& outPrefix);

}  // namespace qe
