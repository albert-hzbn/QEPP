// Harmonic Approximation (HA) thermodynamics computed from a phonon DOS.
//
// All integrals are evaluated with the trapezoidal rule on the frequency
// grid supplied by phonopy (uniform spacing, frequency in THz, cyclic).
//
// Key identities (cyclic frequency ν, not angular ω):
//   x(ν,T) = h·ν / (kB·T)           (dimensionless phonon energy)
//   ZPE     = ½ ∫₀^∞ g(ν)·h·ν dν   [eV/cell]
//   Fvib    = ZPE + kBT ∫ g(ν) ln(1 − e^{−x}) dν
//   S       = kB  ∫ g(ν) [x·n̄ − ln(1 − e^{−x})] dν    (n̄ = Bose–Einstein)
//   Cv      = kB  ∫ g(ν) [x² e^x / (e^x − 1)²] dν
//   Evib    = ZPE + ∫ g(ν)·h·ν·n̄ dν
//
// Physical constants (SI-derived, exact or CODATA-2018):
//   h  = 4.135667696e-3  eV/THz
//   kB = 8.617333262e-5  eV/K
//   NA = 6.02214076e23   mol⁻¹

#include "qe/phonon.hpp"

#include <cmath>
#include <stdexcept>
#include <vector>

namespace qe {

namespace {
constexpr double kPlanck  = 4.135667696e-3;  // h  [eV/THz]
constexpr double kKB      = 8.617333262e-5;  // kB [eV/K]
}  // anonymous namespace

HaResult compute_harmonic_approximation(const PhononDosData& dos,
                                        const std::vector<double>& temperatures,
                                        int natom) {
    const int n = static_cast<int>(dos.freqTHz.size());
    if (n < 2)
        throw std::runtime_error("Phonon DOS has fewer than 2 frequency points.");
    if (dos.totalDos.size() != static_cast<size_t>(n))
        throw std::runtime_error("Phonon DOS: frequency and DOS arrays have different lengths.");

    // Uniform frequency step (trapezoidal — phonopy always outputs uniform grid)
    const double dnu = dos.freqTHz[1] - dos.freqTHz[0];
    if (dnu <= 0.0)
        throw std::runtime_error("Phonon DOS: frequency grid is not monotonically increasing.");

    HaResult result;
    result.natom = natom;

    // ── Zero-point energy ─────────────────────────────────────────────────────
    // ZPE = 0.5 * ∫ g(ν) h ν dν
    // Trapezoidal with endpoint halving:
    double zpe = 0.0;
    for (int i = 0; i < n; ++i) {
        const double nu = dos.freqTHz[i];
        const double g  = dos.totalDos[i];
        if (nu <= 0.0) continue;               // skip ν ≤ 0 (acoustic at Γ)
        const double w = (i == 0 || i == n - 1) ? 0.5 : 1.0;
        zpe += w * g * kPlanck * nu;
    }
    zpe *= 0.5 * dnu;
    result.zeroPointEnergy = zpe;

    // ── Thermal properties ────────────────────────────────────────────────────
    result.thermal.reserve(temperatures.size());

    for (double T : temperatures) {
        HaThermalPoint pt;
        pt.temperature = T;

        // T = 0 K: only ZPE contributes, S = Cv = 0
        if (T < 1e-10) {
            pt.freeEnergy     = zpe;
            pt.entropy        = 0.0;
            pt.heatCapacity   = 0.0;
            pt.internalEnergy = zpe;
            result.thermal.push_back(pt);
            continue;
        }

        const double kT  = kKB * T;
        double fvib = 0.0;   // ∫ g kT ln(1-e^{-x}) dν  (added to ZPE for F)
        double svib = 0.0;   // ∫ g kB [x n̄ − ln(1−e^{−x})] dν
        double cvib = 0.0;   // ∫ g kB x²eˣ/(eˣ−1)² dν
        double evib = 0.0;   // ∫ g hν n̄ dν  (occupation energy, added to ZPE)

        for (int i = 0; i < n; ++i) {
            const double nu = dos.freqTHz[i];
            const double g  = dos.totalDos[i];
            if (nu <= 0.0) continue;           // acoustic Γ: contribution is 0

            const double x   = kPlanck * nu / kT;
            const double w   = (i == 0 || i == n - 1) ? 0.5 : 1.0;

            // Large x: phonon not thermally excited → contributions → 0
            if (x > 500.0) continue;

            const double ex1  = std::expm1(x);           // e^x - 1
            const double bose = 1.0 / ex1;               // n̄ = 1/(e^x - 1)
            const double logf = std::log1p(-std::exp(-x)); // ln(1 - e^{-x})

            fvib += w * g * kT * logf;
            svib += w * g * kKB * (x * bose - logf);
            // Cv = kB * x² * e^x / (e^x - 1)²  =  kB * x² * (bose+1) * bose
            cvib += w * g * kKB * (x * x) * (bose + 1.0) * bose;
            evib += w * g * kPlanck * nu * bose;
        }
        fvib *= dnu;
        svib *= dnu;
        cvib *= dnu;
        evib *= dnu;

        pt.freeEnergy     = zpe + fvib;
        pt.entropy        = svib;
        pt.heatCapacity   = cvib;
        pt.internalEnergy = zpe + evib;
        result.thermal.push_back(pt);
    }

    return result;
}

}  // namespace qe
