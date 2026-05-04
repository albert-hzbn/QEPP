// Temperature-dependent elastic constants via the quasi-harmonic approximation.
//
// Two contributions are summed:
//
//   1. Quasi-static (dominant)
//      C_ij^{QS}(T) = C_ij^{static}(V(T))
//      Obtained by fitting a cubic polynomial to the static elastic constants
//      C_ij^{static}(V_n) computed at N different volumes, then evaluating the
//      polynomial at the QHA equilibrium volume V(T).
//
//   2. Phonon bulk-modulus correction
//      B^{ph}(T) = V(T) · d²F_vib / dV²|_{V(T),T}
//      Computed from a cubic polynomial fit to F_vib(V_n, T) obtained from
//      the phonon DOS at each volume.  Added to the normal–normal components
//      (Voigt rows/cols 0,1,2) of the tensor, keeping shear entries unchanged.
//      This distributes ΔB^{ph} consistently:
//        ΔB_V = (ΔC_11+ΔC_22+ΔC_33+2ΔC_12+2ΔC_13+2ΔC_23)/9
//              = (3 + 6) × B^{ph} / 9 = B^{ph}  ✓

#include "qe/qha_elastic.hpp"
#include "qe/qha.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace qe {

// ── Physical constants ────────────────────────────────────────────────────────
static constexpr double kQhe_RyToEV      = 13.605693122994;
static constexpr double kQhe_KJmolToEV   = 1.0 / 96.485332;   // kJ/mol → eV
static constexpr double kQhe_EVAng3ToGPa = 160.21766208;       // eV/Å³  → GPa

// ── Polynomial utilities ──────────────────────────────────────────────────────

// Fit cubic polynomial y = a0 + a1*x + a2*x² + a3*x³ in a least-squares sense.
static std::array<double, 4> fit_cubic(const std::vector<double>& x,
                                        const std::vector<double>& y) {
    const int n = static_cast<int>(x.size());
    if (n < 4)
        throw std::runtime_error("Cubic fit requires at least 4 points.");
    Eigen::MatrixXd A(n, 4);
    Eigen::VectorXd b(n);
    for (int i = 0; i < n; ++i) {
        A(i, 0) = 1.0;
        A(i, 1) = x[i];
        A(i, 2) = x[i] * x[i];
        A(i, 3) = x[i] * x[i] * x[i];
        b(i)    = y[i];
    }
    const Eigen::VectorXd c = A.colPivHouseholderQr().solve(b);
    return {c[0], c[1], c[2], c[3]};
}

static double eval_poly3(const std::array<double, 4>& c, double x) {
    return c[0] + c[1] * x + c[2] * x * x + c[3] * x * x * x;
}

// Second derivative: d²P/dx² = 2·c[2] + 6·c[3]·x
static double eval_poly3_d2(const std::array<double, 4>& c, double x) {
    return 2.0 * c[2] + 6.0 * c[3] * x;
}

// ── VRH helper ────────────────────────────────────────────────────────────────
static void vrh_from_voigt(const Eigen::Matrix<double, 6, 6>& C,
                            double& KV, double& KR, double& KH,
                            double& GV, double& GR, double& GH,
                            double& EH, double& nuH) {
    Eigen::Matrix<double, 6, 6> S = C.inverse();

    KV = (C(0,0)+C(1,1)+C(2,2) + 2.0*(C(0,1)+C(0,2)+C(1,2))) / 9.0;
    GV = (C(0,0)+C(1,1)+C(2,2) - C(0,1)-C(0,2)-C(1,2)
          + 3.0*(C(3,3)+C(4,4)+C(5,5))) / 15.0;

    const double KRi = S(0,0)+S(1,1)+S(2,2) + 2.0*(S(0,1)+S(0,2)+S(1,2));
    KR = (std::abs(KRi) > 1e-30) ? 1.0/KRi : 0.0;
    const double GRi = 4.0*(S(0,0)+S(1,1)+S(2,2)) - 4.0*(S(0,1)+S(0,2)+S(1,2))
                     + 3.0*(S(3,3)+S(4,4)+S(5,5));
    GR = (std::abs(GRi) > 1e-30) ? 15.0/GRi : 0.0;

    KH  = 0.5 * (KV + KR);
    GH  = 0.5 * (GV + GR);
    EH  = (3.0*KH + GH > 1e-10) ? 9.0*KH*GH / (3.0*KH + GH) : 0.0;
    nuH = (6.0*KH + 2.0*GH > 1e-10) ? (3.0*KH - 2.0*GH) / (2.0*(3.0*KH + GH)) : 0.0;
}

// ── compute_qha_elastic ───────────────────────────────────────────────────────
QhaElasticResult compute_qha_elastic(
    const std::vector<QhaElasticVolumePoint>& elasticVols,
    const std::vector<QhaVolumePoint>&        phonVols) {

    const int N = static_cast<int>(elasticVols.size());
    if (N < 4)
        throw std::runtime_error(
            "QHA elastic requires at least 4 volume points; got " +
            std::to_string(N));
    if (static_cast<int>(phonVols.size()) != N)
        throw std::runtime_error(
            "elastic and phonon volume arrays must have the same length.");

    // ── Collect volumes, static energies ─────────────────────────────────────
    std::vector<double> vols(N), eStaticEV(N);
    for (int i = 0; i < N; ++i) {
        vols[i]      = elasticVols[i].volumeAng3;
        eStaticEV[i] = elasticVols[i].energyRy * kQhe_RyToEV;
    }

    // ── Run standard QHA to get V(T), B_T(T), α(T) ───────────────────────────
    const QhaResult qha = compute_qha(phonVols);

    // ── Fit cubic polynomial for each C_ij(V) element ────────────────────────
    // polyC[r][c] holds the 4 coefficients of C_rc(V).
    std::array<std::array<std::array<double,4>, 6>, 6> polyC;
    for (int r = 0; r < 6; ++r) {
        for (int c = 0; c < 6; ++c) {
            std::vector<double> cij(N);
            for (int n = 0; n < N; ++n)
                cij[n] = elasticVols[n].elastic.C(r, c);
            try {
                polyC[r][c] = fit_cubic(vols, cij);
            } catch (...) {
                polyC[r][c] = {0.0, 0.0, 0.0, 0.0};
            }
        }
    }

    // ── Get crystal family from the first elastic result ──────────────────────
    const std::string family = elasticVols[0].elastic.crystalFamily;

    // ── Assemble result ───────────────────────────────────────────────────────
    QhaElasticResult result;
    result.crystalFamily = family;
    result.staticEos     = qha.staticEos;
    result.nvolumes      = N;

    for (const auto& qhaPt : qha.thermal) {
        const double T   = qhaPt.temperature;
        const double Veq = qhaPt.volumeAng3;

        // ── 1. Quasi-static C at V(T) ─────────────────────────────────────────
        Eigen::Matrix<double, 6, 6> C_qs;
        for (int r = 0; r < 6; ++r)
            for (int c = 0; c < 6; ++c)
                C_qs(r, c) = eval_poly3(polyC[r][c], Veq);

        // Bulk modulus from quasi-static tensor (Voigt)
        const double KV_qs = (C_qs(0,0)+C_qs(1,1)+C_qs(2,2)
                              + 2.0*(C_qs(0,1)+C_qs(0,2)+C_qs(1,2))) / 9.0;

        // ── 2. Phonon bulk-modulus correction ─────────────────────────────────
        // Collect F_vib(V_n, T) [kJ/mol → eV] at this temperature.
        double B_ph = 0.0;
        {
            std::vector<double> fvib(N);
            bool ok = true;
            for (int n = 0; n < N; ++n) {
                const auto it = phonVols[n].fvib.find(T);
                if (it == phonVols[n].fvib.end()) { ok = false; break; }
                fvib[n] = it->second * kQhe_KJmolToEV;  // eV/unit-cell
            }
            if (ok) {
                try {
                    const auto pf = fit_cubic(vols, fvib);
                    // d²F_vib/dV² at Veq [eV/Å⁶], then B^ph = Veq · d²F/dV² [eV/Å³] → GPa
                    const double d2f = eval_poly3_d2(pf, Veq);
                    B_ph = Veq * d2f * kQhe_EVAng3ToGPa;
                } catch (...) {
                    B_ph = 0.0;
                }
            }
        }

        // ── 3. Total C: quasi-static + phonon correction ─────────────────────
        // The isotropic phonon contribution (bulk only, shear unchanged)
        // adds B_ph to all normal–normal Voigt entries (r,c ∈ {0,1,2}).
        // This preserves ΔB_V = (3+6)·B_ph / 9 = B_ph and ΔG = 0.
        Eigen::Matrix<double, 6, 6> C_tot = C_qs;
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c)
                C_tot(r, c) += B_ph;

        // ── 4. VRH moduli from C_total ────────────────────────────────────────
        QhaElasticThermalPoint pt;
        pt.temperature   = T;
        pt.volumeAng3    = Veq;
        pt.alpha         = qhaPt.alpha;
        pt.bulkModStatic = KV_qs;          // Voigt bulk from quasi-static tensor
        pt.bulkModPhonon = B_ph;
        pt.bulkMod       = KV_qs + B_ph;  // total (approximate)
        pt.C_static      = C_qs;
        pt.C_total       = C_tot;

        vrh_from_voigt(C_tot,
                       pt.KV, pt.KR, pt.KH,
                       pt.GV, pt.GR, pt.GH,
                       pt.EH, pt.nuH);

        result.thermal.push_back(pt);
    }

    return result;
}

}  // namespace qe
