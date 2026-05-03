#include "qe/dos.hpp"

#include <matplot/matplot.h>

#include <algorithm>
#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>

namespace qe {

void write_dos_plot_bundle(const std::string& dosInputPath,
                           const std::vector<std::vector<double>>& dosRows,
                           double fermiEv,
                           const std::string& outPrefix) {
    const std::string dataPath = outPrefix + ".dos.dat";
    const std::string pngPath = outPrefix + ".dos.png";

    std::ofstream dataOut(dataPath);
    if (!dataOut.is_open()) {
        throw std::runtime_error("Could not open DOS output data file: " + dataPath);
    }

    dataOut << "# E_minus_Ef(eV)";
    for (size_t c = 1; c < dosRows.front().size(); ++c) {
        dataOut << " col" << (c + 1);
    }
    dataOut << "\n";

    dataOut << std::fixed << std::setprecision(10);
    for (const auto& row : dosRows) {
        dataOut << (row[0] - fermiEv);
        for (size_t c = 1; c < row.size(); ++c) {
            dataOut << " " << row[c];
        }
        dataOut << "\n";
    }
    dataOut.close();

    std::vector<double> x;
    x.reserve(dosRows.size());
    for (const auto& row : dosRows) {
        x.push_back(row[0] - fermiEv);
    }

    using namespace matplot;
    auto fig = figure(true);
    fig->size(1200, 800);
    hold(on);

    const std::vector<std::array<float, 3>> palette = {
        {0.12F, 0.47F, 0.71F}, {0.84F, 0.15F, 0.16F}, {0.17F, 0.63F, 0.17F},
        {0.58F, 0.40F, 0.74F}, {1.00F, 0.50F, 0.05F}};

    const size_t nCols = dosRows.front().size();
    double yMin = std::numeric_limits<double>::infinity();
    double yMax = -std::numeric_limits<double>::infinity();
    for (size_t c = 1; c < nCols; ++c) {
        std::vector<double> y;
        y.reserve(dosRows.size());
        for (const auto& row : dosRows) {
            y.push_back(row[c]);
            yMin = std::min(yMin, row[c]);
            yMax = std::max(yMax, row[c]);
        }

        auto p = plot(x, y);
        p->display_name("DOS col" + std::to_string(c + 1));
        p->line_width(1.8);
        p->color(palette[(c - 1) % palette.size()]);
    }

    auto efLine =
        plot(std::vector<double>{0.0, 0.0}, std::vector<double>{yMin, yMax}, "--k");
    efLine->display_name("E_F");
    title("Quantum ESPRESSO DOS");
    xlabel("E - E_F (eV)");
    ylabel("DOS (states/eV)");
    grid(on);
    legend();

    std::cout << "Prepared DOS data: " << dataPath << "\n";
    save(pngPath);
    std::cout << "Generated DOS plot image (Matplot++): " << pngPath << "\n";
    std::cout << "Source DOS file: " << dosInputPath << "\n";
}

}  // namespace qe

namespace qe {

PdosSummary write_pdos_plots(const std::string& dosPath,
                             const std::vector<std::vector<double>>& totalDosRows,
                             double fermiEv,
                             const std::string& outPrefix) {
    if (totalDosRows.empty()) return {};

    const auto pdos = parse_pdos_summary(dosPath);
    if (pdos.channelCount == 0) {
        std::cout << "  No PDOS files found — looked for '<prefix>.pdos_atm#*'"
                  << " alongside the DOS input file.\n"
                  << "  Run projwfc.x first to generate partial-DOS files.\n";
        return pdos;
    }

    // Total DOS energy axis (for background grey line)
    const size_t nTot = totalDosRows.size();
    std::vector<double> xTot(nTot), tdos(nTot);
    for (size_t i = 0; i < nTot; ++i) {
        xTot[i] = totalDosRows[i][0] - fermiEv;
        tdos[i] = totalDosRows[i][1];
    }

    const size_t nPdos = pdos.energies.size();
    std::vector<double> xPdos(nPdos);
    for (size_t i = 0; i < nPdos; ++i) xPdos[i] = pdos.energies[i] - fermiEv;

    const double yMax = *std::max_element(tdos.begin(), tdos.end());

    const std::vector<std::array<float, 3>> palette = {
        {0.12f, 0.47f, 0.71f}, {0.84f, 0.15f, 0.16f}, {0.17f, 0.63f, 0.17f},
        {0.58f, 0.40f, 0.74f}, {1.00f, 0.50f, 0.05f}, {0.55f, 0.34f, 0.29f},
        {0.89f, 0.47f, 0.76f}, {0.50f, 0.50f, 0.50f}};

    using namespace matplot;

    // ── Plot 1: elemental PDOS ──────────────────────────────────────────────
    {
        const std::string pngPath = outPrefix + ".pdos_elem.png";
        auto fig = figure(true);
        fig->size(1200, 800);
        hold(on);

        // Total DOS as light grey background (may span wider range than PDOS)
        auto pt = plot(xTot, tdos);
        pt->display_name("Total DOS");
        pt->line_width(1.5f);
        pt->color({0.75f, 0.75f, 0.75f});

        size_t ci = 0;
        for (const auto& [elem, v] : pdos.byElement) {
            auto p = plot(xPdos, v);
            p->display_name(elem);
            p->line_width(2.2f);
            p->color(palette[ci++ % palette.size()]);
        }

        auto ef = plot(std::vector<double>{0.0, 0.0},
                       std::vector<double>{0.0, yMax * 1.05}, "--k");
        ef->display_name("E_F");

        title("Elemental PDOS");
        xlabel("E - E_F (eV)");
        ylabel("PDOS (states/eV)");
        grid(on);
        legend();
        save(pngPath);
        std::cout << "Generated elemental PDOS plot : " << pngPath << "\n";
    }

    // ── Plot 2: orbital-type PDOS ───────────────────────────────────────────
    {
        const std::string pngPath = outPrefix + ".pdos_orb.png";
        const std::map<char, std::string> orbLabel = {
            {'s', "s"}, {'p', "p"}, {'d', "d"}, {'f', "f"}};

        auto fig = figure(true);
        fig->size(1200, 800);
        hold(on);

        auto pt = plot(xTot, tdos);
        pt->display_name("Total DOS");
        pt->line_width(1.5f);
        pt->color({0.75f, 0.75f, 0.75f});

        size_t ci = 0;
        for (const auto& [orb, v] : pdos.byOrbital) {
            const std::string lbl = orbLabel.count(orb) ? orbLabel.at(orb)
                                                        : std::string(1, orb);
            auto p = plot(xPdos, v);
            p->display_name(lbl);
            p->line_width(2.2f);
            p->color(palette[ci++ % palette.size()]);
        }

        auto ef = plot(std::vector<double>{0.0, 0.0},
                       std::vector<double>{0.0, yMax * 1.05}, "--k");
        ef->display_name("E_F");

        title("Orbital PDOS");
        xlabel("E - E_F (eV)");
        ylabel("PDOS (states/eV)");
        grid(on);
        legend();
        save(pngPath);
        std::cout << "Generated orbital PDOS plot   : " << pngPath << "\n";
    }

    // ── Summary ─────────────────────────────────────────────────────────────
    std::cout << "PDOS summary: " << pdos.channelCount << " channels  |  "
              << pdos.byElement.size() << " element(s): ";
    for (const auto& [e, _] : pdos.byElement) std::cout << e << " ";
    std::cout << " |  " << pdos.byOrbital.size() << " orbital type(s): ";
    for (const auto& [o, _] : pdos.byOrbital) std::cout << o << " ";
    std::cout << "\n";

    return pdos;
}

}  // namespace qe
