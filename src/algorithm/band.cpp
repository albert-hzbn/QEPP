#include "qe/band.hpp"

#include <algorithm>
#include <iostream>

#include "qe/utils.hpp"

namespace qe {

static std::string orbital_name(int l) {
    switch (l) {
        case 0: return "s";
        case 1: return "p";
        case 2: return "d";
        case 3: return "f";
        default: return "l" + std::to_string(l);
    }
}

std::vector<FatBandGroup> build_fatband_groups(
    const AtomicProj& proj,
    const std::vector<std::string>& elements,
    const std::vector<int>& atomNums,
    const std::vector<int>& orbitalLs) {

    std::vector<FatBandGroup> groups;

    if (!elements.empty() && !orbitalLs.empty()) {
        // One group per (element, orbital) pair
        for (const auto& elem : elements) {
            for (int l : orbitalLs) {
                FatBandGroup g;
                g.name = elem + "-" + orbital_name(l);
                for (int iw = 0; iw < proj.nwfc; ++iw) {
                    const auto& wi = proj.wfcInfo[static_cast<size_t>(iw)];
                    if (to_lower(wi.elem) == to_lower(elem) && wi.l == l)
                        g.wfcIdx.push_back(iw);
                }
                if (!g.wfcIdx.empty()) groups.push_back(std::move(g));
            }
        }
    } else if (!elements.empty()) {
        // One group per element
        for (const auto& elem : elements) {
            FatBandGroup g; g.name = elem;
            for (int iw = 0; iw < proj.nwfc; ++iw) {
                const auto& wi = proj.wfcInfo[static_cast<size_t>(iw)];
                if (to_lower(wi.elem) == to_lower(elem))
                    g.wfcIdx.push_back(iw);
            }
            if (!g.wfcIdx.empty()) groups.push_back(std::move(g));
        }
    } else if (!atomNums.empty() && !orbitalLs.empty()) {
        // One group per (atom, orbital) pair
        for (int anum : atomNums) {
            for (int l : orbitalLs) {
                FatBandGroup g;
                std::string elemName;
                for (int iw = 0; iw < proj.nwfc; ++iw) {
                    const auto& wi = proj.wfcInfo[static_cast<size_t>(iw)];
                    if (wi.atomnum == anum && wi.l == l) {
                        if (elemName.empty()) elemName = wi.elem;
                        g.wfcIdx.push_back(iw);
                    }
                }
                g.name = (elemName.empty() ? "atom" : elemName) + "#"
                         + std::to_string(anum) + "-" + orbital_name(l);
                if (!g.wfcIdx.empty()) groups.push_back(std::move(g));
            }
        }
    } else if (!atomNums.empty()) {
        // One group per atom
        for (int anum : atomNums) {
            FatBandGroup g;
            std::string elemName;
            for (int iw = 0; iw < proj.nwfc; ++iw) {
                const auto& wi = proj.wfcInfo[static_cast<size_t>(iw)];
                if (wi.atomnum == anum) {
                    if (elemName.empty()) elemName = wi.elem;
                    g.wfcIdx.push_back(iw);
                }
            }
            g.name = (elemName.empty() ? "atom" : elemName) + "#" + std::to_string(anum);
            if (!g.wfcIdx.empty()) groups.push_back(std::move(g));
        }
    } else if (!orbitalLs.empty()) {
        // One group per orbital, all atoms summed
        for (int l : orbitalLs) {
            FatBandGroup g;
            g.name = orbital_name(l);
            for (int iw = 0; iw < proj.nwfc; ++iw)
                if (proj.wfcInfo[static_cast<size_t>(iw)].l == l)
                    g.wfcIdx.push_back(iw);
            if (!g.wfcIdx.empty()) groups.push_back(std::move(g));
        }
    } else {
        // Auto: one group per distinct element
        std::vector<std::string> seen;
        for (const auto& wi : proj.wfcInfo) {
            if (std::find(seen.begin(), seen.end(), wi.elem) == seen.end())
                seen.push_back(wi.elem);
        }
        for (const auto& elem : seen) {
            FatBandGroup g; g.name = elem;
            for (int iw = 0; iw < proj.nwfc; ++iw)
                if (proj.wfcInfo[static_cast<size_t>(iw)].elem == elem)
                    g.wfcIdx.push_back(iw);
            if (!g.wfcIdx.empty()) groups.push_back(std::move(g));
        }
    }

    if (groups.empty()) {
        std::cerr << "Warning: no matching atoms/elements found in atomic_proj.xml; "
                     "plotting all contributions.\n";
        FatBandGroup g; g.name = "total";
        for (int iw = 0; iw < proj.nwfc; ++iw) g.wfcIdx.push_back(iw);
        groups.push_back(std::move(g));
    }

    return groups;
}

}  // namespace qe
