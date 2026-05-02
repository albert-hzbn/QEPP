#include "qe/charge.hpp"

#include <matplot/matplot.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>

#include "qe/utils.hpp"

namespace fs = std::filesystem;

namespace qe {

// ============================================================================
// Cube parser
// ============================================================================

CubeData parse_cube(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open()) {
        throw std::runtime_error("Could not open cube file: " + path);
    }

    CubeData cube;
    std::string line;

    // Lines 1-2: comments
    std::getline(in, cube.comment1);
    std::getline(in, cube.comment2);

    // Line 3: natoms  xo  yo  zo
    {
        std::getline(in, line);
        std::istringstream iss(line);
        iss >> cube.natoms >> cube.origin[0] >> cube.origin[1] >> cube.origin[2];
        // natoms < 0 means density is in MO format — handle absolute value
        cube.natoms = std::abs(cube.natoms);
    }

    // Lines 4-6: ni  vx  vy  vz  for each axis
    for (int axis = 0; axis < 3; ++axis) {
        std::getline(in, line);
        std::istringstream iss(line);
        int n = 0;
        double vx = 0, vy = 0, vz = 0;
        iss >> n >> vx >> vy >> vz;
        cube.dims[axis] = n;
        cube.axes(axis, 0) = vx;
        cube.axes(axis, 1) = vy;
        cube.axes(axis, 2) = vz;
    }

    // Lines 7..6+natoms: Z  charge  x  y  z
    cube.atomZ.resize(static_cast<size_t>(cube.natoms));
    cube.atomPos.resize(static_cast<size_t>(cube.natoms));
    for (int i = 0; i < cube.natoms; ++i) {
        std::getline(in, line);
        std::istringstream iss(line);
        double charge = 0;
        iss >> cube.atomZ[i] >> charge
            >> cube.atomPos[i][0] >> cube.atomPos[i][1] >> cube.atomPos[i][2];
    }

    // Volumetric data: NX * NY * NZ values, outer→inner = x, y, z
    const size_t total = static_cast<size_t>(cube.dims[0]) *
                         static_cast<size_t>(cube.dims[1]) *
                         static_cast<size_t>(cube.dims[2]);
    cube.values.reserve(total);
    double v = 0;
    while (in >> v) {
        cube.values.push_back(v);
    }

    if (cube.values.size() != total) {
        throw std::runtime_error(
            "Cube file '" + path + "': expected " + std::to_string(total) +
            " data points but got " + std::to_string(cube.values.size()) + ".");
    }

    return cube;
}

// ============================================================================
// Pre-processing helpers
// ============================================================================

static void write_one_pp_input(const std::string& filePath,
                               const std::string& prefix,
                               const std::string& outdir,
                               int plotNum,
                               const std::string& filplot,
                               const std::string& cubeOut) {
    std::ofstream out(filePath);
    if (!out.is_open()) {
        throw std::runtime_error("Could not write pp.x input: " + filePath);
    }

    out << "&INPUTPP\n"
        << "  prefix     = '" << prefix  << "',\n"
        << "  outdir     = '" << outdir  << "',\n"
        << "  plot_num   = " << plotNum << ",\n"
        << "  filplot    = '" << filplot << "'\n"
        << "/\n"
        << "&PLOT\n"
        << "  nfile          = 1,\n"
        << "  filepp(1)      = '" << filplot << "',\n"
        << "  weight(1)      = 1.0,\n"
        << "  iflag          = 3,\n"
        << "  output_format  = 6,\n"
        << "  fileout        = '" << cubeOut << "'\n"
        << "/\n";
}

void write_pp_inputs(const std::string& scfInputPath,
                     const std::string& outputDir) {
    const auto lines = load_lines(scfInputPath);
    std::string prefix = extract_quoted_assignment(lines, "prefix");
    std::string outdir = extract_quoted_assignment(lines, "outdir");

    if (prefix.empty()) prefix = "qe";
    if (outdir.empty()) outdir = "./tmp";

    fs::create_directories(outputDir);

    const std::string pfx = outputDir + "/" + prefix;

    // Charge density (plot_num = 0)
    write_one_pp_input(pfx + ".charge.pp.in",  prefix, outdir,
                       0, prefix + ".charge",      prefix + ".charge.cube");

    // Charge density difference (plot_num = 6)
    write_one_pp_input(pfx + ".charge_diff.pp.in", prefix, outdir,
                       6, prefix + ".charge_diff",  prefix + ".charge_diff.cube");

    // Electron localisation function (plot_num = 8)
    write_one_pp_input(pfx + ".elf.pp.in",     prefix, outdir,
                       8, prefix + ".elf",          prefix + ".elf.cube");

    std::cout << "Written pp.x inputs in: " << outputDir << "/\n"
              << "  " << prefix << ".charge.pp.in      (plot_num=0  — charge density)\n"
              << "  " << prefix << ".charge_diff.pp.in (plot_num=6  — charge-density difference)\n"
              << "  " << prefix << ".elf.pp.in         (plot_num=8  — ELF)\n"
              << "\n"
              << "Run each with:\n"
              << "  cd " << outputDir << " && pp.x < " << prefix << ".charge.pp.in > "
              << prefix << ".charge.pp.out\n"
              << "  cd " << outputDir << " && pp.x < " << prefix << ".charge_diff.pp.in > "
              << prefix << ".charge_diff.pp.out\n"
              << "  cd " << outputDir << " && pp.x < " << prefix << ".elf.pp.in > "
              << prefix << ".elf.pp.out\n"
              << "\nThen plot with:\n"
              << "  qepp charge -post " << prefix << ".charge.cube\n"
              << "  qepp charge -post " << prefix << ".charge_diff.cube charge_diff\n"
              << "  qepp charge -post " << prefix << ".elf.cube elf\n";
}

// ============================================================================
// Plot helpers
// ============================================================================

// Extract a 2D slice perpendicular to `axis` at integer index `idx`.
// Returns a (n1 x n2) matrix (row-major, matplot convention).
static std::vector<std::vector<double>>
extract_slice(const CubeData& cube, int axis, int idx) {
    const int nx = cube.dims[0];
    const int ny = cube.dims[1];
    const int nz = cube.dims[2];

    if (axis == 2) {
        // XY plane at given z index → rows=y, cols=x
        std::vector<std::vector<double>> mat(static_cast<size_t>(ny),
                                             std::vector<double>(static_cast<size_t>(nx)));
        for (int ix = 0; ix < nx; ++ix)
            for (int iy = 0; iy < ny; ++iy)
                mat[static_cast<size_t>(iy)][static_cast<size_t>(ix)] =
                    cube.at(ix, iy, idx);
        return mat;
    } else if (axis == 1) {
        // XZ plane at given y index → rows=z, cols=x
        std::vector<std::vector<double>> mat(static_cast<size_t>(nz),
                                             std::vector<double>(static_cast<size_t>(nx)));
        for (int ix = 0; ix < nx; ++ix)
            for (int iz = 0; iz < nz; ++iz)
                mat[static_cast<size_t>(iz)][static_cast<size_t>(ix)] =
                    cube.at(ix, idx, iz);
        return mat;
    } else {
        // YZ plane at given x index → rows=z, cols=y
        std::vector<std::vector<double>> mat(static_cast<size_t>(nz),
                                             std::vector<double>(static_cast<size_t>(ny)));
        for (int iy = 0; iy < ny; ++iy)
            for (int iz = 0; iz < nz; ++iz)
                mat[static_cast<size_t>(iz)][static_cast<size_t>(iy)] =
                    cube.at(idx, iy, iz);
        return mat;
    }
}

// For charge_diff: find max absolute value for symmetric colorbar.
static double abs_max(const std::vector<std::vector<double>>& mat) {
    double vmax = 0.0;
    for (const auto& row : mat)
        for (double v : row)
            vmax = std::max(vmax, std::abs(v));
    return vmax;
}

void write_charge_plots(const std::string& cubePath,
                        const std::string& outPrefix,
                        const std::string& quantity) {
    const CubeData cube = parse_cube(cubePath);
    const int nx = cube.dims[0];
    const int ny = cube.dims[1];
    const int nz = cube.dims[2];

    // Determine axis labels and color settings from quantity type
    const bool isDiff = (quantity == "charge_diff");
    const bool isElf  = (quantity == "elf");

    // ── 1. Three orthogonal slices — one square figure per slice ──────────
    {
        using namespace matplot;

        const char* axisLabel[3]  = {"YZ (x=mid)", "XZ (y=mid)", "XY (z=mid)"};
        const char* sliceSuffix[3] = {".slice_YZ.png", ".slice_XZ.png", ".slice_XY.png"};
        const int midIdx[3]        = {nx / 2, ny / 2, nz / 2};

        // extract_slice(cube, panel, idx) [n_row x n_col]:
        //   panel 0 (YZ): col = cube axis 1 (y), row = cube axis 2 (z)
        //   panel 1 (XZ): col = cube axis 0 (x), row = cube axis 2 (z)
        //   panel 2 (XY): col = cube axis 0 (x), row = cube axis 1 (y)
        const int colAxis[3] = {1, 0, 0};
        const int rowAxis[3] = {2, 2, 1};
        const double bohr2ang = 0.529177;

        for (int panel = 0; panel < 3; ++panel) {
            const auto slice = extract_slice(cube, panel, midIdx[panel]);

            // Physical extents in Å
            const double xmax = cube.dims[colAxis[panel]] *
                                cube.axes.row(colAxis[panel]).norm() * bohr2ang;
            const double ymax = cube.dims[rowAxis[panel]] *
                                cube.axes.row(rowAxis[panel]).norm() * bohr2ang;

            // Size the figure so the plot area is physically square:
            // pixel_w / pixel_h = xmax / ymax (plus fixed colorbar margin)
            const int colorbarPx = 120;
            const int basePx     = 600;
            const int figW = static_cast<int>(basePx * (xmax / ymax)) + colorbarPx;
            const int figH = basePx;

            auto fig = figure(true);
            fig->size(static_cast<size_t>(figW), static_cast<size_t>(figH));

            imagesc(0.0, xmax, 0.0, ymax, slice);

            if (isDiff) {
                const double vmax = abs_max(slice);
                if (vmax > 0.0) gca()->zlim({-vmax, vmax});
                colormap(palette::rdbu());
            } else if (isElf) {
                gca()->zlim({0.0, 1.0});
                colormap(palette::viridis());
            } else {
                colormap(palette::plasma());
            }

            colorbar();
            ::matplot::title(axisLabel[panel]);
            xlabel("x (\\AA)");
            ylabel("y (\\AA)");

            const std::string pngPath = outPrefix + sliceSuffix[panel];
            save(pngPath);
            std::cout << "Generated slice plot: " << pngPath << "\n";
        }
    }

    // ── 2. 3-D scatter (isosurface proxy) ──────────────────────────────────
    {
        using namespace matplot;
        auto fig = figure(true);
        fig->size(1000, 900);

        // Compute a threshold: use mean + 1*stddev (or 0.2 for ELF)
        const double total = static_cast<double>(cube.values.size());
        const double mean  = std::accumulate(cube.values.begin(), cube.values.end(), 0.0) / total;

        double threshold = 0.0;
        if (isElf) {
            threshold = 0.75;  // ELF > 0.75 = highly localised electrons
        } else if (isDiff) {
            // Accumulation region: positive lobe
            double posSum  = 0.0;
            int    posCnt  = 0;
            for (double v : cube.values) if (v > 0) { posSum += v; ++posCnt; }
            threshold = (posCnt > 0) ? posSum / posCnt * 0.5 : 0.0;
        } else {
            double var = 0.0;
            for (double v : cube.values) var += (v - mean) * (v - mean);
            var /= total;
            threshold = mean + std::sqrt(var);
        }

        // Convert voxel indices to Angstrom coordinates
        const double bohr2ang = 0.529177;
        std::vector<double> xs, ys, zs, cs;
        xs.reserve(2000); ys.reserve(2000); zs.reserve(2000); cs.reserve(2000);

        // Subsample to keep at most ~8000 scatter points for responsiveness
        const int stride = std::max(1, static_cast<int>(
            std::cbrt(static_cast<double>(cube.values.size()) / 8000.0)));

        for (int ix = 0; ix < nx; ix += stride) {
            for (int iy = 0; iy < ny; iy += stride) {
                for (int iz = 0; iz < nz; iz += stride) {
                    const double v = cube.at(ix, iy, iz);
                    if (v < threshold) continue;

                    // Physical position = origin + ix*ax + iy*ay + iz*az
                    const Eigen::Vector3d pos =
                        (cube.origin +
                         ix * cube.axes.row(0).transpose() +
                         iy * cube.axes.row(1).transpose() +
                         iz * cube.axes.row(2).transpose()) * bohr2ang;

                    xs.push_back(pos[0]);
                    ys.push_back(pos[1]);
                    zs.push_back(pos[2]);
                    cs.push_back(v);
                }
            }
        }

        if (!xs.empty()) {
            scatter3(xs, ys, zs, std::vector<double>(xs.size(), 4.0), cs);
        }

        // Overlay atom positions as larger markers
        {
            std::vector<double> ax_pos, ay_pos, az_pos;
            for (const auto& p : cube.atomPos) {
                ax_pos.push_back(p[0] * bohr2ang);
                ay_pos.push_back(p[1] * bohr2ang);
                az_pos.push_back(p[2] * bohr2ang);
            }
            if (!ax_pos.empty()) {
                hold(on);
                scatter3(ax_pos, ay_pos, az_pos,
                         std::vector<double>(ax_pos.size(), 60.0));
            }
        }

        std::string qtLabel;
        if (isElf)       qtLabel = "ELF > 0.75";
        else if (isDiff) qtLabel = "charge density difference (positive lobe)";
        else             qtLabel = "charge density (above mean + 1σ)";

        title("3-D: " + qtLabel);
        xlabel("x (A)"); ylabel("y (A)"); zlabel("z (A)");
        colormap(isDiff ? palette::hot() : (isElf ? palette::viridis() : palette::plasma()));
        colorbar();

        const std::string scatterPng = outPrefix + ".3d.png";
        save(scatterPng);
        std::cout << "Generated 3-D scatter plot: " << scatterPng << "\n";
    }
}

}  // namespace qe
