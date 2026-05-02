#include "qe/conv.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

#include <matplot/matplot.h>

namespace qe {

void write_conv_report(const std::vector<ConvPoint>& pts,
                       const std::string& paramType,
                       const std::string& outPrefix) {
    if (pts.empty())
        throw std::runtime_error("No convergence data to report.");

    const std::string paramLabel =
        (paramType == "kspacing") ? "k-spacing (1/Ang)" : "ecutwfc (Ry)";

    constexpr int w1 = 20, w2 = 20, w3 = 16;

    auto printTable = [&](std::ostream& os) {
        os << std::left
           << std::setw(w1) << paramLabel
           << std::setw(w2) << "E_tot (eV)"
           << std::setw(w3) << "|dE| (meV/at)\n"
           << std::string(w1 + w2 + w3, '-') << "\n";
        for (const auto& pt : pts) {
            os << std::left
               << std::setw(w1) << std::fixed << std::setprecision(4) << pt.param
               << std::setw(w2) << std::fixed << std::setprecision(6) << pt.energyEv;
            if (pt.deltaEv < 1e-9)
                os << std::setw(w3) << "(ref)";
            else
                os << std::setw(w3) << std::fixed << std::setprecision(3) << pt.deltaEv;
            os << "\n";
        }
    };

    std::cout << "\nConvergence test: " << paramType << "\n";
    printTable(std::cout);

    const std::string txtPath = outPrefix + ".conv.txt";
    std::ofstream ftxt(txtPath);
    if (!ftxt.is_open())
        throw std::runtime_error("Could not create: " + txtPath);
    ftxt << "Convergence test: " << paramType << "\n";
    printTable(ftxt);
    ftxt.close();
    std::cout << "Saved: " << txtPath << "\n";

    // ── PNG: |dE| vs parameter ────────────────────────────────────────────────
    std::vector<double> xs, ys;
    for (const auto& pt : pts) {
        xs.push_back(pt.param);
        ys.push_back(pt.deltaEv);
    }

    using namespace matplot;
    auto fig = figure(true);
    fig->size(700, 450);
    auto ax = fig->current_axes();
    plot(ax, xs, ys, "-o")
        ->color({0.12f, 0.47f, 0.71f})
        .line_width(1.5f)
        .marker_size(6.0f);
    xlabel(ax, paramLabel);
    ylabel(ax, "|\\DeltaE| (meV/atom)");
    title(ax, "Convergence: " + paramType);
    grid(ax, true);

    const std::string pngPath = outPrefix + ".conv.png";
    save(fig, pngPath);
    std::cout << "Saved: " << pngPath << "\n";
}

}  // namespace qe
