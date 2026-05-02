#include "qe/stm.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include "qe/utils.hpp"

namespace fs = std::filesystem;

namespace qe {

void write_stm_pp_input(const std::string& scfInputPath,
                        const std::string& outputDir,
                        double biasEv) {
    const auto lines = load_lines(scfInputPath);
    std::string prefix = extract_quoted_assignment(lines, "prefix");
    std::string outdir = extract_quoted_assignment(lines, "outdir");
    if (prefix.empty()) prefix = "qe";
    if (outdir.empty()) outdir = "./tmp";

    fs::create_directories(outputDir);

    const double emin = (biasEv >= 0.0) ? 0.0 : biasEv;
    const double emax = (biasEv >= 0.0) ? biasEv : 0.0;

    const std::string filplot  = prefix + ".stm";
    const std::string cubeOut  = prefix + ".stm.cube";
    const std::string ppInPath = outputDir + "/" + prefix + ".stm.pp.in";

    std::ofstream out(ppInPath);
    if (!out.is_open()) throw std::runtime_error("Could not create: " + ppInPath);

    out << "&INPUTPP\n"
        << "  prefix     = '" << prefix << "',\n"
        << "  outdir     = '" << outdir << "',\n"
        << "  plot_num   = 5,\n"
        << "  emin       = " << emin << ",\n"
        << "  emax       = " << emax << ",\n"
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

    std::cout << "Written:  " << ppInPath << "\n";
    std::cout << "Bias:     [" << emin << ", " << emax << "] eV relative to E_F\n";
    std::cout << "Run:      pp.x < " << ppInPath << " > " << prefix << ".stm.pp.out\n";
    std::cout << "Then:     qepp stm -post " << cubeOut << "\n";
}

}  // namespace qe
