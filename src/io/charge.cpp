#include "qe/charge.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "qe/utils.hpp"

namespace fs = std::filesystem;

namespace qe {

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
        iss >> cube.atomZ[i] >> charge >> cube.atomPos[i][0] >> cube.atomPos[i][1] >> cube.atomPos[i][2];
    }

    // Volumetric data: NX * NY * NZ values, outer→inner = x, y, z
    const size_t total = static_cast<size_t>(cube.dims[0]) * static_cast<size_t>(cube.dims[1]) *
                         static_cast<size_t>(cube.dims[2]);
    cube.values.reserve(total);
    double v = 0;
    while (in >> v) {
        cube.values.push_back(v);
    }

    if (cube.values.size() != total) {
        throw std::runtime_error("Cube file '" + path + "': expected " + std::to_string(total) +
                                 " data points but got " + std::to_string(cube.values.size()) + ".");
    }

    return cube;
}

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
        << "  prefix     = '" << prefix << "',\n"
        << "  outdir     = '" << outdir << "',\n"
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

    if (prefix.empty()) {
        prefix = "qe";
    }
    if (outdir.empty()) {
        outdir = "./tmp";
    }

    fs::create_directories(outputDir);

    const std::string pfx = outputDir + "/" + prefix;

    // Charge density (plot_num = 0)
    write_one_pp_input(pfx + ".charge.pp.in", prefix, outdir, 0, prefix + ".charge", prefix + ".charge.cube");

    // Charge density difference (plot_num = 6)
    write_one_pp_input(pfx + ".charge_diff.pp.in", prefix, outdir, 6, prefix + ".charge_diff",
                       prefix + ".charge_diff.cube");

    // Electron localisation function (plot_num = 8)
    write_one_pp_input(pfx + ".elf.pp.in", prefix, outdir, 8, prefix + ".elf", prefix + ".elf.cube");

    std::cout << "Written pp.x inputs in: " << outputDir << "/\n"
              << "  " << prefix << ".charge.pp.in      (plot_num=0  — charge density)\n"
              << "  " << prefix << ".charge_diff.pp.in (plot_num=6  — charge-density difference)\n"
              << "  " << prefix << ".elf.pp.in         (plot_num=8  — ELF)\n"
              << "\n"
              << "Run each with:\n"
              << "  cd " << outputDir << " && pp.x < " << prefix << ".charge.pp.in > " << prefix
              << ".charge.pp.out\n"
              << "  cd " << outputDir << " && pp.x < " << prefix << ".charge_diff.pp.in > " << prefix
              << ".charge_diff.pp.out\n"
              << "  cd " << outputDir << " && pp.x < " << prefix << ".elf.pp.in > " << prefix << ".elf.pp.out\n"
              << "\nThen plot with:\n"
              << "  qepp charge -post " << prefix << ".charge.cube\n"
              << "  qepp charge -post " << prefix << ".charge_diff.cube charge_diff\n"
              << "  qepp charge -post " << prefix << ".elf.cube elf\n";
}

}  // namespace qe
