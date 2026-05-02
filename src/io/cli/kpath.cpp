#include "qe/cli/kpath.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "qe/cif.hpp"
#include "qe/help.hpp"

namespace qe {

int handle_kpath_mode(int argc, char** argv, int s) {
    if (argc < 3 + s || argc > 5 + s) {
        print_help_command(argv[0], "kpath", "-pre");
        return 1;
    }

    const std::string cifPath = argv[2 + s];
    const int pointsPerSegment = (argc >= 4 + s) ? std::stoi(argv[3 + s]) : 20;

    const CifStructure structure = parse_cif(cifPath);
    const SymmetryKPath kpath = suggest_kpath_from_cif(structure);

    std::cout << "Crystal family  : " << kpath.family << "\n";
    if (!structure.spaceGroupName.empty() || structure.spaceGroupNumber > 0) {
        std::cout << "Space group     : "
                  << (structure.spaceGroupName.empty() ? "?" : structure.spaceGroupName)
                  << "  (IT# " << structure.spaceGroupNumber << ")\n";
    }
    std::cout << "K-path          :";
    for (size_t i = 0; i < kpath.nodes.size(); ++i) {
        std::cout << " " << kpath.nodes[i].label;
        if (i + 1 < kpath.nodes.size()) std::cout << " -";
    }
    std::cout << "\n\n";

    std::ostringstream block;
    block << "K_POINTS crystal_b\n";
    block << "  " << kpath.nodes.size() << "\n";
    block << std::fixed << std::setprecision(8);
    for (size_t i = 0; i < kpath.nodes.size(); ++i) {
        const int w = (i + 1 < kpath.nodes.size()) ? pointsPerSegment : 1;
        block << "  " << std::setw(11) << kpath.nodes[i].k.x()
              << " " << std::setw(11) << kpath.nodes[i].k.y()
              << " " << std::setw(11) << kpath.nodes[i].k.z()
              << " " << w << " ! " << kpath.nodes[i].label << "\n";
    }
    std::cout << block.str();

    if (argc >= 5 + s) {
        const std::string outPath = argv[4 + s];
        std::ofstream out(outPath);
        if (!out.is_open()) {
            throw std::runtime_error("Cannot open output file: " + outPath);
        }
        out << block.str();
        std::cout << "\nWritten to: " << outPath << "\n";
    }

    return 0;
}

}  // namespace qe
