#include "qe/cli/commands.hpp"

#include <Eigen/Dense>

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "qe/cif.hpp"
#include "qe/help.hpp"
#include "qe/qe_input.hpp"
#include "qe/types.hpp"
#include "qe/utils.hpp"

namespace qe {

int handle_cif_mode(int argc, char** argv, int s) {
    if (argc < 4 + s || argc > 7 + s) {
        print_help_command(argv[0], "cif", "-pre");
        return 1;
    }

    const std::string cifPath = argv[2 + s];
    const double kspacing = std::stod(argv[3 + s]);
    if (kspacing <= 0.0) {
        throw std::runtime_error("kspacing must be > 0.");
    }

    const std::string defaultOutput = stem_from_path(cifPath) + ".scf.in";
    const std::string outPath = (argc >= 5 + s) ? argv[4 + s] : defaultOutput;
    const int ecutwfc = (argc >= 6 + s) ? std::stoi(argv[5 + s]) : 50;
    const int ecutrho = (argc >= 7 + s) ? std::stoi(argv[6 + s]) : (8 * ecutwfc);

    CifStructure structure = parse_cif(cifPath);
    std::vector<std::string> species = build_species_blocks(structure.atoms);
    const Eigen::Vector3i kGrid = kmesh_from_kspacing(structure.cellAngstrom, kspacing);
    const Eigen::Vector3i kShift(0, 0, 0);

    const std::string prefix = stem_from_path(cifPath);
    write_qe_input(outPath, prefix, structure, ecutwfc, ecutrho, kGrid, kShift,
                   species, kspacing);

    std::cout << "Generated Quantum ESPRESSO input file: " << outPath << "\n";
    std::cout << "Detected " << structure.atoms.size() << " atoms, " << species.size()
              << " species\n";
    std::cout << "Generated k-mesh: " << kGrid.x() << " " << kGrid.y() << " "
              << kGrid.z() << "\n";
    if (!structure.spaceGroupName.empty() || structure.spaceGroupNumber > 0) {
        std::cout << "Symmetry: "
                  << (structure.spaceGroupName.empty() ? "(name not found)"
                                                       : structure.spaceGroupName)
                  << "  IT#: " << structure.spaceGroupNumber << "\n";
    }
    return 0;
}

}  // namespace qe
