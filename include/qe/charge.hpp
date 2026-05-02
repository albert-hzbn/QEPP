#pragma once

#include <Eigen/Dense>

#include <array>
#include <string>
#include <vector>

namespace qe {

// ---------------------------------------------------------------------------
// Gaussian cube file data (used by QE pp.x with output_format=6)
// ---------------------------------------------------------------------------
struct CubeData {
    std::string comment1;
    std::string comment2;

    int natoms = 0;
    Eigen::Vector3d origin{0.0, 0.0, 0.0};  // Bohr
    std::array<int, 3> dims{};               // {nx, ny, nz}
    Eigen::Matrix3d axes;   // row i = step vector for axis i, in Bohr

    std::vector<int>           atomZ;
    std::vector<Eigen::Vector3d> atomPos;  // Bohr

    // values[ix * ny * nz + iy * nz + iz]  (outer→inner: x, y, z)
    std::vector<double> values;

    double at(int ix, int iy, int iz) const {
        return values[static_cast<size_t>(ix) * static_cast<size_t>(dims[1]) *
                          static_cast<size_t>(dims[2]) +
                      static_cast<size_t>(iy) * static_cast<size_t>(dims[2]) +
                      static_cast<size_t>(iz)];
    }
};

// ---------------------------------------------------------------------------
// Parse a Gaussian cube file (.cube) produced by QE pp.x
// ---------------------------------------------------------------------------
CubeData parse_cube(const std::string& path);

// ---------------------------------------------------------------------------
// Pre-processing: write three pp.x input files for
//   charge density (plot_num=0), charge-density difference (plot_num=6),
//   and electron localisation function (plot_num=8).
// scfInputPath is used to extract prefix, outdir, pseudo_dir from the SCF run.
// outputDir is where the pp.x inputs are written (created if absent).
// ---------------------------------------------------------------------------
void write_pp_inputs(const std::string& scfInputPath,
                     const std::string& outputDir);

// ---------------------------------------------------------------------------
// Post-processing: read a cube file and produce plots.
//   outPrefix   — output filename prefix  (e.g. "si_charge")
//   quantity    — label used in titles / colour-bar labels
//                 "charge"      → plasma colormap, no diverging scale
//                 "charge_diff" → RdBu colormap, symmetric around zero
//                 "elf"         → viridis colormap, clamped to [0, 1]
// Outputs:
//   <outPrefix>.slice.png   — three orthogonal midplane heatmaps (x, y, z)
//   <outPrefix>.3d.png      — 3-D scatter of high-value points (isosurface proxy)
// ---------------------------------------------------------------------------
void write_charge_plots(const std::string& cubePath,
                        const std::string& outPrefix,
                        const std::string& quantity = "charge");

}  // namespace qe
