#pragma once

#include <Eigen/Dense>

#include <string>
#include <vector>

#include "qe/types.hpp"

namespace qe {

double atomic_mass(const std::string& symbol);
Eigen::Vector3i kmesh_from_kspacing(const Eigen::Matrix3d& cellAngstrom,
                                    double kspacingInvA);
std::vector<std::string> build_species_blocks(const std::vector<Atom>& atoms);

void write_qe_input(const std::string& fileName,
                    const std::string& prefix,
                    const CifStructure& structure,
                    int ecutwfc,
                    int ecutrho,
                    const Eigen::Vector3i& kGrid,
                    const Eigen::Vector3i& kShift,
                    const std::vector<std::string>& speciesBlocks,
                    double kspacingInvA);

void write_bands_input_from_scf_template(const std::string& scfInputPath,
                                         const std::string& bandsInputPath,
                                         const SymmetryKPath& kpath,
                                         int pointsPerSegment,
                                         int nbnd = 0);

void write_bands_pp_input_from_scf_template(const std::string& scfInputPath,
                                            const std::string& bandsPpInputPath,
                                            const std::string& filbandName);

}  // namespace qe
