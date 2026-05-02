#pragma once

#include <string>

#include "qe/types.hpp"

namespace qe {

// ── Pre-processing ────────────────────────────────────────────────────────────
// Generate deformed QE SCF inputs for the energy-strain method.
// Creates:
//   <outDir>/
//     <pattern>/
//       d<+/-><delta>/
//         <prefix>.in
//
// scfTemplatePath : existing QE SCF input (ibrav=0, CELL_PARAMETERS angstrom)
// outDir          : root output directory (created if missing)
// nDeltas         : number of strain points per pattern (odd, ≥ 5)
// maxDelta        : maximum strain magnitude (e.g. 0.04)
void generate_elastic_inputs(const std::string& scfTemplatePath,
                              const std::string& outDir,
                              int nDeltas = 7,
                              double maxDelta = 0.04);

// ── Post-processing ───────────────────────────────────────────────────────────
// Collect total energies from QE output files under <outDir>, fit E(δ) curves,
// extract the full 6×6 stiffness matrix (using crystal symmetry detected from
// the template), compute derived properties and print a report.
//
// scfTemplatePath : same template used during pre-processing (cell + mass)
// outDir          : same root used during pre-processing
// Returns         : filled ElasticResults struct
ElasticResults compute_elastic_properties(const std::string& scfTemplatePath,
                                          const std::string& outDir);

// Print a formatted report to stdout and optionally to a text file.
void print_elastic_report(const ElasticResults& res,
                          const std::string& crystalFamily,
                          const std::string& outFile = "");

}  // namespace qe
