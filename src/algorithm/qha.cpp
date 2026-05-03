// Quasi-Harmonic Approximation numerical algorithms.
//
// References:
//   Phonopy QHA: https://phonopy.github.io/phonopy/qha.html
//   Birch-Murnaghan EOS: F. Birch, Phys. Rev. 71, 809 (1947)
//   Thermal expansion: A. Togo, I. Tanaka, Scr. Mater. 108, 1-5 (2015)
//
// Physical constants used throughout:
//   1 Ry  = 13.605693122994 eV
//   1 eV  = 96.485332 kJ/mol
//   GPa conversion: 1 eV/Ang^3 = 160.21766208 GPa

#include "qe/qha.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace qe {

static constexpr double kRyToEV    = 13.605693122994;
static constexpr double kEVToKJmol = 96.485332;
static constexpr double kKJmolToEV = 1.0 / kEVToKJmol;
static constexpr double kJKmolToEV = kKJmolToEV / 1000.0;  // J/K/mol → eV/K
static constexpr double kEVAng3ToGPa = 160.21766208;       // eV/Ang^3 → GPa
static constexpr double kNAvogadro = 6.02214076e23;
static constexpr double kKB        = 8.617333262e-5;        // eV/K (Boltzmann)

// ── Birch-Murnaghan 3rd-order EOS ─────────────────────────────────────────────
// E(V) = E0 + (9*V0*B0/16) * { [(V0/V)^(2/3)-1]^3 * B0'
//              + [(V0/V)^(2/3)-1]^2 * [6-4*(V0/V)^(2/3)] }
// All energies in eV, volumes in Ang^3, B in GPa (converted internally).

double bm_eos_energy(double V, const BmEosParams& p) {
    const double x = std::pow(p.V0 / V, 2.0 / 3.0) - 1.0;
    const double y = std::pow(p.V0 / V, 2.0 / 3.0);
    // 9*V0*B0/16, B0 in GPa → eV/Ang^3
    const double coeff = 9.0 * p.V0 * (p.B0 / kEVAng3ToGPa) / 16.0;
    return p.E0 + coeff * (x * x * x * p.B0p + x * x * (6.0 - 4.0 * y));
}

double bm_eos_pressure(double V, const BmEosParams& p) {
    // P = (3*B0/2)*[(V0/V)^(7/3)-(V0/V)^(5/3)]*{1+(3/4)*(B0'-4)*[(V0/V)^(2/3)-1]}
    const double r = std::pow(p.V0 / V, 1.0 / 3.0);  // (V0/V)^(1/3)
    const double r7 = r * r * r * r * r * r * r;      // (V0/V)^(7/3)
    const double r5 = r * r * r * r * r;               // (V0/V)^(5/3)
    const double x  = r * r - 1.0;                    // (V0/V)^(2/3) - 1
    const double P_eV_A3 = 1.5 * (p.B0 / kEVAng3ToGPa) * (r7 - r5) *
                           (1.0 + 0.75 * (p.B0p - 4.0) * x);
    return P_eV_A3 * kEVAng3ToGPa;  // GPa
}

double bm_eos_bulk_modulus(double V, const BmEosParams& p) {
    // B_T = V * d^2E/dV^2 (in eV/Ang^3), converted to GPa
    // Use finite difference for simplicity
    const double dV = V * 1e-5;
    const double Ep = bm_eos_energy(V + dV, p);
    const double Em = bm_eos_energy(V - dV, p);
    const double E0 = bm_eos_energy(V, p);
    const double d2EdV2 = (Ep - 2.0 * E0 + Em) / (dV * dV);
    return V * d2EdV2 * kEVAng3ToGPa;
}

// ── Nelder-Mead simplex minimiser (4-parameter, no external deps) ─────────────
// Minimises f(x) where x is a 4-vector.

using Vec4 = std::array<double, 4>;

static double nm_eval(const Vec4& x,
                      const std::vector<double>& vols,
                      const std::vector<double>& energies) {
    BmEosParams p;
    p.E0  = x[0];
    p.V0  = x[1];
    p.B0  = std::abs(x[2]);   // must be positive
    p.B0p = x[3];
    double rms = 0.0;
    for (size_t i = 0; i < vols.size(); ++i) {
        const double diff = bm_eos_energy(vols[i], p) - energies[i];
        rms += diff * diff;
    }
    return rms;
}

static Vec4 nm_centroid(const std::vector<Vec4>& simplex, int worst) {
    Vec4 c = {0, 0, 0, 0};
    int n = static_cast<int>(simplex.size());
    for (int i = 0; i < n; ++i) {
        if (i == worst) continue;
        for (int k = 0; k < 4; ++k) c[k] += simplex[i][k];
    }
    for (int k = 0; k < 4; ++k) c[k] /= static_cast<double>(n - 1);
    return c;
}

static Vec4 nm_reflect(const Vec4& c, const Vec4& worst, double alpha) {
    Vec4 r;
    for (int k = 0; k < 4; ++k) r[k] = c[k] + alpha * (c[k] - worst[k]);
    return r;
}

BmEosParams fit_bm_eos(const std::vector<double>& volumes,
                        const std::vector<double>& energies) {
    if (volumes.size() < 4)
        throw std::runtime_error("BM EOS fit requires at least 4 (V, E) points.");
    if (volumes.size() != energies.size())
        throw std::runtime_error("Volume and energy arrays must have the same length.");

    // Initial guess: E0 = min E, V0 = V at min E, B0 = 100 GPa, B0' = 4
    size_t imin = static_cast<size_t>(
        std::min_element(energies.begin(), energies.end()) - energies.begin());
    double E0_init = energies[imin];
    double V0_init = volumes[imin];

    // Rough B0 from curvature around minimum: B ≈ V * d²E/dV²
    double B0_init = 100.0;
    if (imin > 0 && imin < volumes.size() - 1) {
        const double dV1 = volumes[imin + 1] - volumes[imin];
        const double dV0 = volumes[imin] - volumes[imin - 1];
        if (std::abs(dV1) > 1e-10 && std::abs(dV0) > 1e-10) {
            const double d2E = 2.0 * ((energies[imin + 1] - energies[imin]) / dV1 -
                                      (energies[imin] - energies[imin - 1]) / dV0) /
                               (dV1 + dV0);
            const double Braw = V0_init * d2E * kEVAng3ToGPa;
            if (Braw > 1.0 && Braw < 1000.0) B0_init = Braw;
        }
    }

    // Build 5-vertex simplex (4+1)
    std::vector<Vec4> simplex(5);
    simplex[0] = {E0_init, V0_init, B0_init, 4.0};
    simplex[1] = {E0_init + 0.01, V0_init, B0_init, 4.0};
    simplex[2] = {E0_init, V0_init * 1.02, B0_init, 4.0};
    simplex[3] = {E0_init, V0_init, B0_init * 1.1, 4.0};
    simplex[4] = {E0_init, V0_init, B0_init, 5.0};

    auto f = [&](const Vec4& x) { return nm_eval(x, volumes, energies); };
    std::vector<double> fval(5);
    for (int i = 0; i < 5; ++i) fval[i] = f(simplex[i]);

    const double ftol = 1e-14;
    const int maxIter = 100000;

    for (int iter = 0; iter < maxIter; ++iter) {
        // Sort
        std::vector<int> idx = {0, 1, 2, 3, 4};
        std::sort(idx.begin(), idx.end(), [&](int a, int b) { return fval[a] < fval[b]; });
        const int best  = idx[0];
        const int worst = idx[4];
        const int sw    = idx[3];  // second worst

        if (fval[best] < ftol) break;
        if (iter > 0 && std::abs(fval[worst] - fval[best]) < ftol * (1.0 + std::abs(fval[best]))) break;

        Vec4 cen = nm_centroid(simplex, worst);

        // Reflection
        Vec4 xr = nm_reflect(cen, simplex[worst], 1.0);
        double fr = f(xr);

        if (fr < fval[best]) {
            // Expansion
            Vec4 xe = nm_reflect(cen, simplex[worst], 2.0);
            double fe = f(xe);
            if (fe < fr) { simplex[worst] = xe; fval[worst] = fe; }
            else          { simplex[worst] = xr; fval[worst] = fr; }
        } else if (fr < fval[sw]) {
            simplex[worst] = xr;
            fval[worst] = fr;
        } else {
            // Contraction
            Vec4 xc = nm_reflect(cen, simplex[worst], 0.5);
            double fc = f(xc);
            if (fc < fval[worst]) {
                simplex[worst] = xc;
                fval[worst] = fc;
            } else {
                // Shrink
                for (int i = 1; i < 5; ++i) {
                    for (int k = 0; k < 4; ++k)
                        simplex[idx[i]][k] = simplex[best][k] +
                                              0.5 * (simplex[idx[i]][k] - simplex[best][k]);
                    fval[idx[i]] = f(simplex[idx[i]]);
                }
            }
        }
    }

    // Pick best vertex
    int bi = static_cast<int>(
        std::min_element(fval.begin(), fval.end()) - fval.begin());
    const Vec4& best = simplex[bi];

    BmEosParams result;
    result.E0  = best[0];
    result.V0  = best[1];
    result.B0  = std::abs(best[2]);
    result.B0p = best[3];
    result.rms = std::sqrt(fval[bi] / static_cast<double>(volumes.size()));
    return result;
}

// ── QHA core computation ──────────────────────────────────────────────────────
// Algorithm (follows phonopy-qha methodology):
//
//  1. For each T, build F(V) = E_static(V) + F_vib(V,T)
//  2. Fit F(V) with BM EOS → V_eq(T), G(T) = F(V_eq), B_T(T)
//  3. Thermal expansion α = d ln V / dT  (central difference in T)
//  4. Cp = Cv + T*V*α²*B_T / N_A  (Nernst-Lindemann correction)
//  5. Grüneisen γ = V*α*B_T / Cv

QhaResult compute_qha(const std::vector<QhaVolumePoint>& points) {
    if (points.size() < 4)
        throw std::runtime_error(
            "QHA requires at least 4 volume points; got " +
            std::to_string(points.size()));

    // ── Collect static energies and fit T=0 EOS ───────────────────────────
    const int N = static_cast<int>(points.size());
    const int natom = points[0].natom;

    std::vector<double> vols(N), eStaticEV(N);
    for (int i = 0; i < N; ++i) {
        vols[i]       = points[i].volumeAng3;
        eStaticEV[i]  = points[i].energyRy * kRyToEV;
    }

    QhaResult result;
    result.natom         = natom;
    result.nvolumes      = N;
    result.inputVolumes  = vols;
    result.inputEnergies = eStaticEV;
    result.staticEos     = fit_bm_eos(vols, eStaticEV);

    // ── Find common temperature grid ──────────────────────────────────────
    // Intersect temperature keys across all volume points.
    std::vector<double> temps;
    {
        // Collect all T values from the first point, then keep only those
        // present in every other point.
        if (points[0].fvib.empty())
            throw std::runtime_error("No phonon free energy data in volume point 0.");

        for (const auto& [T, _] : points[0].fvib)
            temps.push_back(T);

        for (int i = 1; i < N; ++i) {
            std::vector<double> keep;
            for (double T : temps) {
                if (points[i].fvib.count(T)) keep.push_back(T);
            }
            temps = keep;
        }
    }
    if (temps.empty())
        throw std::runtime_error("No common temperature points across all volume files.");
    std::sort(temps.begin(), temps.end());

    // ── Compute QHA at each temperature ──────────────────────────────────
    // We store V_eq(T) to compute α numerically later.
    struct ThermalRaw {
        double T, Veq, G, F, BT, Cv, S;
    };
    std::vector<ThermalRaw> raw;
    raw.reserve(temps.size());

    for (double T : temps) {
        // F(V) = E_static(V) + F_vib(V, T)  [all in eV]
        std::vector<double> ftotal(N);
        double cvAvg = 0.0, sAvg = 0.0;
        for (int i = 0; i < N; ++i) {
            const double fvib_kJ  = points[i].fvib.at(T);
            ftotal[i] = eStaticEV[i] + fvib_kJ * kKJmolToEV;  // eV/unit-cell
            // Use volume-averaged Cv and S at each T (weighted by proximity to V_eq)
            if (points[i].cvib.count(T)) cvAvg += points[i].cvib.at(T);
            if (points[i].svib.count(T)) sAvg  += points[i].svib.at(T);
        }
        cvAvg /= static_cast<double>(N);
        sAvg  /= static_cast<double>(N);

        // Fit BM EOS to F(V,T)
        BmEosParams eos;
        try {
            eos = fit_bm_eos(vols, ftotal);
        } catch (...) {
            continue;  // Skip temperatures where fit fails
        }

        ThermalRaw tr;
        tr.T   = T;
        tr.Veq = eos.V0;                                        // Ang^3
        tr.F   = eos.E0;                                        // eV/unit-cell (Helmholtz at V_eq)
        tr.G   = eos.E0;                                        // at P=0: G = F (Gibbs ≈ Helmholtz)
        tr.BT  = eos.B0;                                        // GPa
        // Use the phonopy Cv/S values (interpolated at V_eq, here averaged)
        tr.Cv  = cvAvg / static_cast<double>(natom);            // J/K/mol/atom
        tr.S   = sAvg  / static_cast<double>(natom);            // J/K/mol/atom
        raw.push_back(tr);
    }

    if (raw.size() < 3)
        throw std::runtime_error("Too few temperature points survived QHA fitting.");

    // ── Compute thermal expansion α = d(ln V)/dT (central differences) ───
    const int M = static_cast<int>(raw.size());
    for (int j = 0; j < M; ++j) {
        double alpha = 0.0;
        if (j == 0) {
            if (M > 1) {
                const double dT = raw[j + 1].T - raw[j].T;
                if (std::abs(dT) > 1e-10)
                    alpha = (raw[j + 1].Veq - raw[j].Veq) / (raw[j].Veq * dT);
            }
        } else if (j == M - 1) {
            const double dT = raw[j].T - raw[j - 1].T;
            if (std::abs(dT) > 1e-10)
                alpha = (raw[j].Veq - raw[j - 1].Veq) / (raw[j].Veq * dT);
        } else {
            const double dT = raw[j + 1].T - raw[j - 1].T;
            if (std::abs(dT) > 1e-10)
                alpha = (raw[j + 1].Veq - raw[j - 1].Veq) / (raw[j].Veq * dT);
        }

        QhaThermalPoint pt;
        pt.temperature       = raw[j].T;
        pt.volumeAng3        = raw[j].Veq;
        pt.volumeAng3PerAtom = raw[j].Veq / static_cast<double>(natom);
        pt.helmholtz         = raw[j].F / static_cast<double>(natom);   // eV/atom
        pt.gibbs             = raw[j].G / static_cast<double>(natom);   // eV/atom
        pt.bulkMod           = raw[j].BT;
        pt.alpha             = alpha;
        pt.cv                = raw[j].Cv;
        pt.entropy           = raw[j].S;

        // Cp correction: Cp = Cv + T * V * α² * B_T / N_A
        // Cv in J/K/mol/atom, V in m^3/atom, B in Pa
        // V_atom = V_ang3_per_atom * 1e-30 m^3
        // B_Pa = B_GPa * 1e9
        const double V_m3  = pt.volumeAng3PerAtom * 1.0e-30;
        const double B_Pa  = pt.bulkMod * 1.0e9;
        const double cp_correction_JKmol =
            raw[j].T * V_m3 * alpha * alpha * B_Pa * kNAvogadro;
        pt.cp = pt.cv + cp_correction_JKmol;

        // Grüneisen parameter γ = V * α * B_T / Cv
        // Cv in J/K/m^3: Cv_m3 = Cv_JKmol/atom / V_m3 / N_A
        const double cv_per_m3 = (pt.cv > 1e-12) ?
            (pt.cv / kNAvogadro / V_m3) : 0.0;
        pt.gruneisen = (cv_per_m3 > 1e-30) ?
            (alpha * B_Pa / cv_per_m3) : 0.0;

        result.thermal.push_back(pt);
    }

    return result;
}

}  // namespace qe
