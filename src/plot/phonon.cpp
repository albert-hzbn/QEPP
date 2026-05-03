// Phonon DOS and band structure plots using Matplot++.

#include "qe/phonon.hpp"

#include <matplot/matplot.h>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace qe {

// ── Color palette ─────────────────────────────────────────────────────────────
static const std::vector<std::array<float, 3>> kPalette = {
    {0.84F, 0.15F, 0.16F},  // red
    {0.12F, 0.47F, 0.71F},  // blue
    {0.17F, 0.63F, 0.17F},  // green
    {0.58F, 0.40F, 0.74F},  // purple
    {1.00F, 0.50F, 0.05F},  // orange
    {0.09F, 0.75F, 0.81F},  // cyan
    {0.85F, 0.80F, 0.10F},  // yellow-green
    {0.50F, 0.30F, 0.10F},  // brown
};

// ── Phonon DOS plot ───────────────────────────────────────────────────────────

void write_phonon_dos_plot(const PhononDosData& dos, const std::string& outPrefix) {
    if (dos.freqTHz.empty())
        throw std::runtime_error("Phonon DOS data is empty.");

    // Save data file
    const std::string datPath = outPrefix + ".phonon_dos.dat";
    const std::string pngPath = outPrefix + ".phonon_dos.png";

    {
        std::ofstream f(datPath);
        if (!f.is_open()) throw std::runtime_error("Cannot write: " + datPath);
        f << "# Phonon DOS\n# freq(THz)   total_dos";
        for (const auto& lbl : dos.projLabels) f << "   pdos_" << lbl;
        f << "\n";
        for (size_t i = 0; i < dos.freqTHz.size(); ++i) {
            f << std::fixed << std::setprecision(8)
              << dos.freqTHz[i] << "   " << dos.totalDos[i];
            for (size_t a = 0; a < dos.projDos.size(); ++a)
                if (i < dos.projDos[a].size()) f << "   " << dos.projDos[a][i];
            f << "\n";
        }
    }
    std::cout << "Saved: " << datPath << "\n";

    // Summary info
    {
        const double fmin = *std::min_element(dos.freqTHz.begin(), dos.freqTHz.end());
        const double fmax = *std::max_element(dos.freqTHz.begin(), dos.freqTHz.end());
        std::cout << std::fixed << std::setprecision(3);
        std::cout << "  Frequency range : " << fmin << " — " << fmax << " THz\n";
        std::cout << "  Frequency points: " << dos.freqTHz.size() << "\n";
        if (dos.hasProjected)
            std::cout << "  Projected atoms : " << dos.projLabels.size() << "\n";
    }

    // ── Matplot++ ────────────────────────────────────────────────────────────
    using namespace matplot;
    auto fig = figure(true);
    fig->size(1000, 700);
    hold(on);

    const double dosMax = *std::max_element(dos.totalDos.begin(), dos.totalDos.end());
    const double fMin   = dos.freqTHz.front();
    const double fMax   = dos.freqTHz.back();

    // Zero-frequency line (to mark acoustic branches)
    if (fMin < 0.0) {
        auto vl = plot(std::vector<double>{0.0, 0.0},
                       std::vector<double>{0.0, dosMax * 1.05}, "--");
        vl->color({0.5F, 0.5F, 0.5F, 1.0F});
        vl->line_width(1.0);
    }

    // Projected DOS (stacked, lighter)
    if (dos.hasProjected) {
        for (size_t a = 0; a < dos.projDos.size(); ++a) {
            const auto& pd = dos.projDos[a];
            const size_t n = std::min(pd.size(), dos.freqTHz.size());
            std::vector<double> fx(dos.freqTHz.begin(), dos.freqTHz.begin() + n);
            std::vector<double> fy(pd.begin(), pd.begin() + n);
            auto p = plot(fx, fy);
            p->display_name(dos.projLabels[a]);
            p->line_width(1.4);
            const auto& c = kPalette[a % kPalette.size()];
            p->color({c[0], c[1], c[2], 0.75F});
        }
    }

    // Total DOS (solid, on top)
    {
        auto p = plot(dos.freqTHz, dos.totalDos);
        p->display_name("Total DOS");
        p->line_width(2.0);
        p->color({0.1F, 0.1F, 0.1F, 1.0F});
    }

    xlim({fMin, fMax});
    ylim({0.0, dosMax * 1.08});

    title("Phonon Density of States");
    xlabel("Frequency (THz)");
    ylabel("DOS (states/THz)");
    if (dos.hasProjected) legend();
    grid(on);

    save(pngPath);
    std::cout << "Saved: " << pngPath << "\n";
}

// ── Phonon band structure plot ────────────────────────────────────────────────

void write_phonon_band_plot(const PhononBandData& band, const std::string& outPrefix) {
    if (band.distances.empty() || band.freqByBand.empty())
        throw std::runtime_error("Phonon band data is empty.");

    const std::string datPath = outPrefix + ".phonon_band.dat";
    const std::string pngPath = outPrefix + ".phonon_band.png";

    // Save data file
    {
        std::ofstream f(datPath);
        if (!f.is_open()) throw std::runtime_error("Cannot write: " + datPath);
        f << "# Phonon band structure\n# dist   freq_band1(THz)  freq_band2 ...\n";
        for (int q = 0; q < band.nqpoints; ++q) {
            f << std::fixed << std::setprecision(8) << band.distances[q];
            for (int b = 0; b < band.nbands; ++b)
                f << "   " << band.freqByBand[b][q];
            f << "\n";
        }
    }
    std::cout << "Saved: " << datPath << "\n";

    // Summary info
    {
        std::cout << "  Q-points : " << band.nqpoints << "\n";
        std::cout << "  Bands    : " << band.nbands   << "\n";
        if (!band.labelMarks.empty()) {
            std::cout << "  Labels   :";
            for (const auto& [d, lbl] : band.labelMarks) std::cout << "  " << lbl;
            std::cout << "\n";
        }
        double fmin = 0.0;
        for (int b = 0; b < band.nbands; ++b)
            for (double f : band.freqByBand[b]) fmin = std::min(fmin, f);
        if (fmin < -0.1)
            std::cout << "  WARNING: imaginary modes present (min freq = "
                      << fmin << " THz)\n";
    }

    // Compute axis limits
    const double kMin = band.distances.front();
    const double kMax = band.distances.back();
    double fMin =  std::numeric_limits<double>::infinity();
    double fMax = -std::numeric_limits<double>::infinity();
    for (int b = 0; b < band.nbands; ++b) {
        for (double f : band.freqByBand[b]) {
            fMin = std::min(fMin, f);
            fMax = std::max(fMax, f);
        }
    }
    const double fPad = (fMax - fMin) * 0.03 + 0.2;
    const double yLo  = fMin - fPad;
    const double yHi  = fMax + fPad;

    // ── Matplot++ ────────────────────────────────────────────────────────────
    using namespace matplot;
    auto fig = figure(true);
    fig->size(1200, 800);
    hold(on);

    // Vertical lines at high-symmetry points
    for (const auto& [d, lbl] : band.labelMarks) {
        if (d > kMin + 1e-10 && d < kMax - 1e-10) {
            auto vl = plot(std::vector<double>{d, d},
                           std::vector<double>{yLo, yHi}, "-");
            vl->color({0.55F, 0.55F, 0.55F, 1.0F});
            vl->line_width(0.8);
        }
    }

    // Zero-frequency dashed line (acoustic)
    auto efLine = plot(std::vector<double>{kMin, kMax},
                       std::vector<double>{0.0, 0.0}, "--");
    efLine->color({0.3F, 0.3F, 0.3F, 0.6F});
    efLine->line_width(1.0);

    // Band curves
    for (int b = 0; b < band.nbands; ++b) {
        auto p = plot(band.distances, band.freqByBand[b]);
        p->line_width(1.6);
        p->color({0.10F, 0.24F, 0.56F, 1.0F});
    }

    // High-symmetry x-tick labels
    if (!band.labelMarks.empty()) {
        std::vector<double>      tickPos;
        std::vector<std::string> tickLbl;
        for (const auto& [d, lbl] : band.labelMarks) {
            tickPos.push_back(d);
            // Convert "GAMMA" → Greek letter representation
            std::string disp = lbl;
            if (disp == "GAMMA" || disp == "Gamma") disp = "{/Symbol G}";
            tickLbl.push_back(disp);
        }
        xticks(tickPos);
        xticklabels(tickLbl);
    }

    xlim({kMin, kMax});
    ylim({yLo, yHi});

    title("Phonon Band Structure");
    xlabel("");
    ylabel("Frequency (THz)");
    grid(off);

    save(pngPath);
    std::cout << "Saved: " << pngPath << "\n";
}

// ── HA thermodynamics plot ────────────────────────────────────────────────────
// Three-panel figure: Cv(T), S(T), F_vib(T).
// Y-axes in J/mol/K  (Cv, S) and eV/atom (F_vib), all per atom.

void write_ha_plot(const HaResult& ha, const std::string& outPrefix) {
    if (ha.thermal.empty())
        throw std::runtime_error("HaResult has no thermal data to plot.");

    const std::string pngPath = outPrefix + ".phonon_ha.png";

    static constexpr double kEVKtoJmolK = 96485.33;  // eV/K/atom → J/mol/K

    const int na = ha.natom;
    std::vector<double> Tvec, Cv_J, S_J, F_eV;

    for (const auto& pt : ha.thermal) {
        if (pt.temperature < 1e-10) continue;
        Tvec.push_back(pt.temperature);
        Cv_J.push_back((pt.heatCapacity  / na) * kEVKtoJmolK);
        S_J .push_back((pt.entropy       / na) * kEVKtoJmolK);
        F_eV.push_back( pt.freeEnergy    / na);
    }

    if (Tvec.empty()) return;
    const double Tmax = Tvec.back();

    // 3R reference line for Dulong-Petit limit (J/mol/K)
    const double threeR = 3.0 * 8.314462618;

    using namespace matplot;
    auto fig = figure(true);
    fig->size(1300, 500);

    // ── Panel 1: Cv ──────────────────────────────────────────────────────────
    subplot(1, 3, 0);
    hold(on);
    {
        auto pCv = plot(Tvec, Cv_J);
        pCv->line_width(2.0);
        pCv->color({0.12F, 0.47F, 0.71F, 1.0F});
        pCv->display_name("Cv (HA)");
    }
    // Dulong-Petit dashed line
    {
        auto dp = plot(std::vector<double>{0.0, Tmax},
                       std::vector<double>{threeR, threeR}, "--");
        dp->color({0.5F, 0.5F, 0.5F, 0.8F});
        dp->line_width(1.0);
        dp->display_name("3R");
    }
    xlim({0.0, Tmax});
    ylim({0.0, threeR * 1.1});
    title("Heat Capacity");
    xlabel("T (K)");
    ylabel("Cv (J/mol/K)");
    grid(on);
    legend();

    // ── Panel 2: S ───────────────────────────────────────────────────────────
    subplot(1, 3, 1);
    hold(on);
    {
        auto pS = plot(Tvec, S_J);
        pS->line_width(2.0);
        pS->color({0.84F, 0.15F, 0.16F, 1.0F});
        pS->display_name("S (HA)");
    }
    xlim({0.0, Tmax});
    title("Entropy");
    xlabel("T (K)");
    ylabel("S (J/mol/K)");
    grid(on);

    // ── Panel 3: F_vib ───────────────────────────────────────────────────────
    subplot(1, 3, 2);
    hold(on);
    {
        auto pF = plot(Tvec, F_eV);
        pF->line_width(2.0);
        pF->color({0.17F, 0.63F, 0.17F, 1.0F});
        pF->display_name("F_vib (HA)");
    }
    xlim({0.0, Tmax});
    title("Vibrational Free Energy");
    xlabel("T (K)");
    ylabel("F_vib (eV/atom)");
    grid(on);

    save(pngPath);
    std::cout << "Saved: " << pngPath << "\n";
}

}  // namespace qe
