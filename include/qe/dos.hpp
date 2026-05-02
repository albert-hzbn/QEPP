#pragma once

#include <string>
#include <vector>

namespace qe {

std::vector<std::vector<double>> parse_dos_table(const std::string& dosPath);
double extract_fermi_from_qe_output(const std::string& qeOutPath);

// Plot total DOS. Always called.
void write_dos_plot_bundle(const std::string& dosInputPath,
                           const std::vector<std::vector<double>>& dosRows,
                           double fermiEv,
                           const std::string& outPrefix);

// Plot elemental and orbital partial DOS from projwfc.x output files.
// Looks for <prefix>.pdos_atm#* files alongside dosPath.
// Silently skips if no PDOS files are present.
void write_pdos_plots(const std::string& dosPath,
                      const std::vector<std::vector<double>>& totalDosRows,
                      double fermiEv,
                      const std::string& outPrefix);

}  // namespace qe
