#include "qe/dos.hpp"

#include <algorithm>
#include <array>
#include <cmath>

namespace qe {
namespace {

struct Moments {
    double m0 = 0.0;
    double m1 = 0.0;
    double m2 = 0.0;
};

void add_trapezoid_moments(Moments& acc,
                           double e0, double d0,
                           double e1, double d1) {
    const double dx = e1 - e0;
    if (dx <= 0.0) return;

    acc.m0 += 0.5 * (d0 + d1) * dx;
    acc.m1 += 0.5 * (e0 * d0 + e1 * d1) * dx;
    acc.m2 += 0.5 * (e0 * e0 * d0 + e1 * e1 * d1) * dx;
}

DBandMethodEstimate estimate_from_moments(const Moments& mom,
                                          double fermiEv,
                                          const std::string& label) {
    DBandMethodEstimate out;
    out.method = label;
    out.integratedWeight = mom.m0;

    constexpr double kTol = 1e-12;
    if (mom.m0 <= kTol) {
        return out;
    }

    const double center = mom.m1 / mom.m0;
    const double variance = std::max(0.0, mom.m2 / mom.m0 - center * center);

    out.valid = true;
    out.centerEv = center;
    out.centerMinusEfEv = center - fermiEv;
    out.widthEv = std::sqrt(variance);
    return out;
}

std::vector<std::array<double, 2>> sorted_nonnegative_ddos(const std::vector<double>& energies,
                                                            const std::vector<double>& dDos) {
    std::vector<std::array<double, 2>> out;
    if (energies.size() != dDos.size()) return out;

    out.reserve(energies.size());
    for (size_t i = 0; i < energies.size(); ++i) {
        out.push_back({energies[i], std::max(0.0, dDos[i])});
    }
    std::sort(out.begin(), out.end(),
              [](const auto& a, const auto& b) { return a[0] < b[0]; });
    return out;
}

}  // namespace

DBandMetrics estimate_d_band_metrics(const PdosSummary& pdos, double fermiEv) {
    DBandMetrics out;
    out.oldMethod.method = "old (occupied states up to E_F)";
    out.newMethod.method = "new (full d-band moment over all states)";

    const auto itD = pdos.byOrbital.find('d');
    if (itD == pdos.byOrbital.end()) {
        return out;
    }

    const auto data = sorted_nonnegative_ddos(pdos.energies, itD->second);
    if (data.size() < 2) {
        return out;
    }
    out.hasDOrbital = true;

    Moments oldMom;
    Moments newMom;

    for (size_t i = 0; i + 1 < data.size(); ++i) {
        const double e0 = data[i][0];
        const double d0 = data[i][1];
        const double e1 = data[i + 1][0];
        const double d1 = data[i + 1][1];

        if (e1 <= e0) {
            continue;
        }

        add_trapezoid_moments(newMom, e0, d0, e1, d1);

        if (e1 <= fermiEv) {
            add_trapezoid_moments(oldMom, e0, d0, e1, d1);
            continue;
        }

        if (e0 < fermiEv && fermiEv < e1) {
            const double t = (fermiEv - e0) / (e1 - e0);
            const double dF = d0 + t * (d1 - d0);
            add_trapezoid_moments(oldMom, e0, d0, fermiEv, dF);
        }
    }

    out.oldMethod = estimate_from_moments(oldMom, fermiEv, out.oldMethod.method);
    out.newMethod = estimate_from_moments(newMom, fermiEv, out.newMethod.method);
    return out;
}

}  // namespace qe
