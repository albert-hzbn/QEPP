#pragma once

#include <map>
#include <string>
#include <vector>

namespace qe {

std::vector<std::vector<double>> parse_dos_table(const std::string& dosPath);
double extract_fermi_from_qe_output(const std::string& qeOutPath);

struct PdosSummary {
    std::vector<double> energies;
    std::map<std::string, std::vector<double>> byElement;
    std::map<char, std::vector<double>> byOrbital;
    size_t channelCount = 0;
};

PdosSummary parse_pdos_summary(const std::string& dosPath);

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
