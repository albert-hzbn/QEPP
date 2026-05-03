#pragma once

// DFPT (Density Functional Perturbation Theory) workflow helpers.
//
// Pre-processing  (phonon -pre):
//   Reads a pw.x SCF input file and writes the four QE input files needed
//   for a full DFPT phonon calculation:
//     1. ph.in       – ph.x phonon calculation (ldisp, q-mesh)
//     2. q2r.in      – q2r.x real-space IFC conversion
//     3. matdyn_dos.in  – matdyn.x phonon DOS
//     4. matdyn_band.in – matdyn.x phonon dispersion (band structure)
//
// Post-processing (phonon -post / phonon -dos / phonon -band / phonon -ha):
//   Reads matdyn.x output files and converts to the same internal
//   PhononDosData / PhononBandData structs used by the phonopy parsers,
//   so all plotting and HA routines work identically for both workflows.
//
// DFPT theory brief:
//   DFPT computes the dynamical matrix D(q) as the second derivative of the
//   DFT total energy with respect to phonon perturbations via the Sternheimer
//   equation.  ph.x diagonalises D(q) at each q-point to obtain phonon
//   frequencies ω(q).  q2r.x Fourier-transforms D(q) to real-space force
//   constants Φ(R), and matdyn.x interpolates Φ(R) back to arbitrary q′ for
//   dispersion and DOS.
//
// Unit convention:
//   matdyn.x writes frequencies in cm⁻¹ by default.  All parse_matdyn_*
//   functions auto-detect the unit from the file header and convert to THz
//   (the same unit used by phonopy) so the rest of the code is oblivious to
//   the source.
//   1 THz (cyclic) = 33.35641 cm⁻¹   [c = 2.99792458e10 cm/s]

#include <string>
#include <vector>

#include "qe/phonon.hpp"  // PhononDosData, PhononBandData

namespace qe {

// ── Pre-processing options ────────────────────────────────────────────────────

struct DfptOptions {
    // q-mesh for ph.x (DFPT phonon calculation)
    int nq1 = 4, nq2 = 4, nq3 = 4;

    // q-mesh for matdyn.x phonon DOS (usually finer than ph.x mesh)
    int nqDos1 = 8, nqDos2 = 8, nqDos3 = 8;

    // ph.x convergence threshold (tr2_ph); 1e-14 gives accurate IFCs
    double tr2_ph = 1.0e-14;

    // Acoustic Sum Rule type for q2r.x and matdyn.x
    // Options: "no", "simple", "crystal", "one-dim", "zero-dim"
    std::string asr = "simple";

    // Whether to include Born effective charges / dielectric constant
    // (epsil = .true. in ph.x; needed for polar insulators, LO-TO splitting)
    bool epsil = false;
};

// Generate the four QE input files needed for a full DFPT phonon calculation.
//
// Reads prefix, outdir, nat, ntyp from scfInputPath.  Writes to outDir/:
//   ph.in, q2r.in, matdyn_dos.in, matdyn_band.in
//
// Returns the paths of all generated files.
std::vector<std::string> generate_phonon_inputs(
    const std::string& scfInputPath,
    const std::string& outDir,
    const DfptOptions& opts = {});

// ── Post-processing parsers ───────────────────────────────────────────────────

// Parse matdyn.x phonon DOS output (fldos file).
//
// Supported formats:
//   - QE default: columns = freq(cm⁻¹)  dos(states/cm⁻¹)  [pdos1  pdos2 ...]
//   - THz mode:   columns = freq(THz)    dos(states/THz)    [pdos1  pdos2 ...]
//
// Units are detected from header comments ("THz" → keep; "cm" or default → convert).
// The returned PhononDosData has frequencies in THz and DOS in states/THz,
// identical to what parse_phonon_dos returns for phonopy total_dos.dat.
PhononDosData parse_matdyn_dos(const std::string& dosPath);

// Parse matdyn.x phonon band output (flfrq file, &plot format).
//
// Format written by matdyn.x with q_in_band_form=.true.:
//   &plot
//     nbnd=N, nks=M /
//   q1x  q1y  q1z
//   f1  f2  ...  fN          (6 per line, in cm⁻¹ or THz)
//   q2x  q2y  q2z
//   ...
//
// qLabelPath: optional file with one line per high-symmetry point:
//   label  qx  qy  qz
// These are matched to q-points in the file by proximity (< 1e-4 tolerance)
// and stored as labelMarks in the returned PhononBandData.
//
// x-axis distances are computed as cumulative |Δq| in crystal coordinates.
// Frequencies are converted to THz if in cm⁻¹.
PhononBandData parse_matdyn_freq(const std::string& freqPath,
                                 const std::string& qLabelPath = "");

// Detect whether a file is matdyn.x DOS output (returns true) or phonopy
// total_dos.dat (returns false).  Heuristic: looks for "cm" or "&plot" in
// the first few lines; falls back to false if uncertain.
bool is_matdyn_dos_format(const std::string& path);

// Detect whether a file is matdyn.x &plot band output (returns true) or a
// phonopy band.yaml file (returns false).
bool is_matdyn_freq_format(const std::string& path);

}  // namespace qe
