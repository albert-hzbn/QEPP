#include "qe/band.hpp"

#include <matplot/matplot.h>

#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace qe {

void write_band_plot_bundle(const std::string& bandInputPath,
                            const BandData& bandData,
                            double fermiEv,
                            const std::string& outPrefix) {
    const std::string dataPath = outPrefix + ".band.dat";
    const std::string pngPath  = outPrefix + ".band.png";

    std::ofstream dataOut(dataPath);
    if (!dataOut.is_open()) {
        throw std::runtime_error("Could not open band output data file: " + dataPath);
    }

    dataOut << std::fixed << std::setprecision(10);
    for (size_t b = 0; b < bandData.kByBand.size(); ++b) {
        const auto& k = bandData.kByBand[b];
        const auto& e = bandData.eByBand[b];
        for (size_t i = 0; i < k.size(); ++i) {
            dataOut << k[i] << " " << (e[i] - fermiEv) << "\n";
        }
        dataOut << "\n";
    }
    dataOut.close();

    using namespace matplot;
    auto fig = figure(true);
    fig->size(1400, 900);
    hold(on);

    double kMin = std::numeric_limits<double>::infinity();
    double kMax = -std::numeric_limits<double>::infinity();
    double eMin = std::numeric_limits<double>::infinity();
    double eMax = -std::numeric_limits<double>::infinity();

    for (size_t b = 0; b < bandData.kByBand.size(); ++b) {
        const auto& k = bandData.kByBand[b];
        for (size_t i = 0; i < k.size(); ++i) {
            kMin = std::min(kMin, k[i]);
            kMax = std::max(kMax, k[i]);
            const double e = bandData.eByBand[b][i] - fermiEv;
            eMin = std::min(eMin, e);
            eMax = std::max(eMax, e);
        }
    }

    const double yPad = (eMax - eMin) * 0.05 + 0.5;
    const double yLo  = std::max(eMin - yPad, -15.0);
    const double yHi  = std::min(eMax + yPad,  15.0);

    if (!bandData.kLabelMarks.empty()) {
        for (const auto& [kPos, label] : bandData.kLabelMarks) {
            if (kPos > kMin + 1e-10 && kPos < kMax - 1e-10) {
                auto vLine = plot(std::vector<double>{kPos, kPos},
                                  std::vector<double>{yLo, yHi}, "--");
                vLine->color({0.55F, 0.55F, 0.55F, 1.0F});
                vLine->line_width(1.0);
            }
        }
    }

    auto efLine = plot(std::vector<double>{kMin, kMax}, std::vector<double>{0.0, 0.0}, "--");
    efLine->color({0.2F, 0.2F, 0.2F, 1.0F});
    efLine->line_width(1.0);

    for (size_t b = 0; b < bandData.kByBand.size(); ++b) {
        const auto& k = bandData.kByBand[b];
        std::vector<double> eShifted;
        eShifted.reserve(bandData.eByBand[b].size());
        for (size_t i = 0; i < k.size(); ++i) {
            eShifted.push_back(bandData.eByBand[b][i] - fermiEv);
        }
        auto p = plot(k, eShifted);
        p->line_width(1.8);
        p->color({0.10F, 0.24F, 0.56F, 1.0F});
    }

    if (!bandData.kLabelMarks.empty()) {
        std::vector<double> tickPositions;
        std::vector<std::string> tickLabels;
        for (const auto& [kPos, label] : bandData.kLabelMarks) {
            tickPositions.push_back(kPos);
            std::string displayLabel = label;
            if (displayLabel == "G" || displayLabel == "GM" || displayLabel == "Gamma") {
                displayLabel = "{/Symbol G}";
            }
            tickLabels.push_back(displayLabel);
        }
        xticks(tickPositions);
        xticklabels(tickLabels);
    }

    xlim({kMin, kMax});
    ylim({yLo, yHi});

    title("Quantum ESPRESSO Band Structure");
    xlabel("");
    ylabel("E - E_{F} (eV)");
    grid(on);

    std::cout << "Prepared band data: " << dataPath << "\n";
    save(pngPath);
    std::cout << "Generated band plot image (Matplot++): " << pngPath << "\n";
    std::cout << "Source band file: " << bandInputPath << "\n";
}

void write_fatband_plots(const BandData& bandData,
                         const std::vector<FatBandGroup>& groups,
                         double fermiEv,
                         const std::string& outPrefix) {
    if (bandData.kByBand.empty()) {
        throw std::runtime_error("Band data is empty.");
    }
    if (groups.empty()) {
        throw std::runtime_error("No projection groups to plot.");
    }

    const int nbnd = static_cast<int>(bandData.kByBand.size());
    const int nk   = static_cast<int>(bandData.kByBand[0].size());

    double kMin =  std::numeric_limits<double>::infinity();
    double kMax = -std::numeric_limits<double>::infinity();
    double eMin =  std::numeric_limits<double>::infinity();
    double eMax = -std::numeric_limits<double>::infinity();
    for (int b = 0; b < nbnd; ++b) {
        for (size_t i = 0; i < bandData.kByBand[b].size(); ++i) {
            kMin = std::min(kMin, bandData.kByBand[b][i]);
            kMax = std::max(kMax, bandData.kByBand[b][i]);
            const double e = bandData.eByBand[b][i] - fermiEv;
            eMin = std::min(eMin, e);
            eMax = std::max(eMax, e);
        }
    }
    const double yPad = (eMax - eMin) * 0.05 + 0.5;
    const double yLo  = std::max(eMin - yPad, -15.0);
    const double yHi  = std::min(eMax + yPad,  15.0);

    // Color palette {r, g, b}
    static const std::array<std::array<float, 3>, 8> kPalette = {{
        {0.85F, 0.15F, 0.15F},
        {0.15F, 0.45F, 0.85F},
        {0.10F, 0.60F, 0.10F},
        {0.85F, 0.55F, 0.05F},
        {0.60F, 0.10F, 0.75F},
        {0.05F, 0.70F, 0.70F},
        {0.85F, 0.80F, 0.10F},
        {0.50F, 0.50F, 0.50F},
    }};

    const double fatLineWidth = 2.0;
    const std::string suffix  = groups.size() == 1 ? groups[0].name : "combined";
    const std::string pngPath = outPrefix + ".fatband_" + suffix + ".png";

    using namespace matplot;
    auto fig = figure(true);
    fig->size(1400, 900);
    hold(on);

    // Vertical lines at high-symmetry points
    for (const auto& [kPos, lbl] : bandData.kLabelMarks) {
        if (kPos > kMin + 1e-10 && kPos < kMax - 1e-10) {
            auto vl = plot(std::vector<double>{kPos, kPos},
                           std::vector<double>{yLo, yHi}, "--");
            vl->color({0.55F, 0.55F, 0.55F, 1.0F});
            vl->line_width(1.0);
        }
    }

    // Fermi level
    auto ef = plot(std::vector<double>{kMin, kMax},
                   std::vector<double>{0.0, 0.0}, "--");
    ef->color({0.3F, 0.3F, 0.3F, 1.0F});
    ef->line_width(1.0);

    // Gray background band lines
    for (int b = 0; b < nbnd; ++b) {
        std::vector<double> es;
        es.reserve(bandData.eByBand[b].size());
        for (double ev : bandData.eByBand[b]) es.push_back(ev - fermiEv);
        auto p = plot(bandData.kByBand[b], es);
        p->color({0.75F, 0.75F, 0.75F, 1.0F});
        p->line_width(1.0);
    }

    // Fat band colored lines — one color per group, uniform line width
    for (size_t gi = 0; gi < groups.size(); ++gi) {
        const auto& grp = groups[gi];
        const auto& col = kPalette[gi % kPalette.size()];
        bool namedGroup = false;

        for (int b = 0; b < nbnd; ++b) {
            const int nkUsed = std::min(nk, static_cast<int>(bandData.kByBand[b].size()));
            std::vector<double> kv, ev;
            kv.reserve(static_cast<size_t>(nkUsed));
            ev.reserve(static_cast<size_t>(nkUsed));
            for (int ik = 0; ik < nkUsed; ++ik) {
                kv.push_back(bandData.kByBand[b][static_cast<size_t>(ik)]);
                ev.push_back(bandData.eByBand[b][static_cast<size_t>(ik)] - fermiEv);
            }
            auto p = plot(kv, ev);
            p->line_width(fatLineWidth);
            p->color(std::array<float, 4>{0.0F, col[0], col[1], col[2]});
            if (!namedGroup) {
                p->display_name(grp.name);
                namedGroup = true;
            } else {
                p->display_name("");
            }
        }
    }

    // X-axis labels at high-symmetry points
    if (!bandData.kLabelMarks.empty()) {
        std::vector<double> tpos;
        std::vector<std::string> tlbl;
        for (const auto& [kPos, lbl] : bandData.kLabelMarks) {
            tpos.push_back(kPos);
            std::string d = lbl;
            if (d == "G" || d == "GM" || d == "Gamma") d = "{/Symbol G}";
            tlbl.push_back(d);
        }
        xticks(tpos);
        xticklabels(tlbl);
    }

    xlim({kMin, kMax});
    ylim({yLo, yHi});
    ylabel("E - E_{F} (eV)");
    grid(on);
    ::matplot::legend();
    gca()->legend()->strings({});

    save(pngPath);
    std::cout << "Generated fat band plot: " << pngPath << "\n";
}

}  // namespace qe
