#include "qe/stm.hpp"

#include <matplot/matplot.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

#include "qe/charge.hpp"

namespace qe {

void write_stm_plot(const std::string& cubePath,
                    const std::string& outPrefix,
                    double heightAng) {
    const CubeData cube = parse_cube(cubePath);
    const int nx = cube.dims[0];
    const int ny = cube.dims[1];
    const int nz = cube.dims[2];

    int zIdx;
    if (heightAng < 0.0) {
        zIdx = nz / 2;
    } else {
        const double dzBohr = cube.axes.row(2).norm();
        zIdx = static_cast<int>(std::round((heightAng / 0.529177) / dzBohr));
        zIdx = std::max(0, std::min(zIdx, nz - 1));
    }

    // Extract XY plane at chosen z-index (rows = y, cols = x)
    std::vector<std::vector<double>> mat(
        static_cast<size_t>(ny), std::vector<double>(static_cast<size_t>(nx)));
    for (int ix = 0; ix < nx; ++ix)
        for (int iy = 0; iy < ny; ++iy)
            mat[static_cast<size_t>(iy)][static_cast<size_t>(ix)] = cube.at(ix, iy, zIdx);

    const double bohr2ang = 0.529177;
    const double xmax = nx * cube.axes.row(0).norm() * bohr2ang;
    const double ymax = ny * cube.axes.row(1).norm() * bohr2ang;
    const double zAng = zIdx * cube.axes.row(2).norm() * bohr2ang;

    using namespace matplot;
    const int colorbarPx = 120;
    const int basePx     = 600;
    const int figW = static_cast<int>(basePx * (xmax / ymax)) + colorbarPx;

    auto fig = figure(true);
    fig->size(static_cast<size_t>(figW), static_cast<size_t>(basePx));

    imagesc(0.0, xmax, 0.0, ymax, mat);
    colormap(palette::hot());
    colorbar();

    const std::string titleStr =
        "STM ILDOS  z = " + std::to_string(zAng).substr(0, 5) + " Ang";
    title(titleStr);
    xlabel("x (Ang)");
    ylabel("y (Ang)");

    const std::string pngPath = outPrefix + ".stm.png";
    save(pngPath);
    std::cout << "Saved: " << pngPath << "  (z = " << zAng << " Ang)\n";
}

}  // namespace qe
