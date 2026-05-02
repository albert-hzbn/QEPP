#include "qe/help.hpp"

#include <iostream>

void print_help(const char* prog) {
    std::cout <<
        "qepp — Quantum ESPRESSO post-processor and input generator\n"
        "\n"
        "USAGE\n"
        "  " << prog << " <command> -pre  [arguments]   (pre-processing)\n"
        "  " << prog << " <command> -post [arguments]   (post-processing)\n"
        "  " << prog << " help [command] [-pre|-post]\n"
        "\n"
        "COMMANDS\n"
        "  cif     -pre   Generate a QE SCF input from a CIF file\n"
        "  dos     -post  Plot density of states from a QE dos.x output\n"
        "  band    -pre   Generate QE band-structure inputs (pw.x + bands.x) from a CIF\n"
        "  band    -post  Plot band structure from a QE bands.x output\n"
        "  kpath   -pre   Suggest a high-symmetry k-path from a CIF file\n"
        "  elastic -pre   Generate deformed SCF inputs for elastic constants\n"
        "  elastic -post  Collect energies and compute elastic constants\n"
        "  charge  -pre   Generate pp.x inputs for charge density / charge diff / ELF\n"
        "  charge  -post  Plot volumetric data from a cube file\n"
        "\n"
        "Run '" << prog << " help <command>' or '" << prog << " <command> --help' for details.\n";
}

void print_help_command(const char* prog, const std::string& cmd,
                        const std::string& sub) {
    // ── cif ──────────────────────────────────────────────────────────────────
    if (cmd == "cif") {
        std::cout <<
            "DESCRIPTION\n"
            "  Reads a CIF (Crystallographic Information File) and generates a ready-to-run\n"
            "  Quantum ESPRESSO pw.x SCF input. The cell vectors, atomic positions, and species\n"
            "  are taken directly from the CIF. A Monkhorst-Pack k-mesh is built automatically\n"
            "  from the requested k-point spacing. Pseudopotential filenames are left as\n"
            "  placeholders that you fill in before running pw.x.\n"
            "\n"
            "USAGE\n"
            "  " << prog << " cif -pre <input.cif> <kspacing> [output.in] [ecutwfc] [ecutrho]\n"
            "\n"
            "ARGUMENTS\n"
            "  input.cif    CIF structure file\n"
            "  kspacing     K-point spacing in 1/Angstrom (e.g. 0.15); smaller = denser mesh\n"
            "  output.in    Output QE input file (default: <stem>.scf.in)\n"
            "  ecutwfc      Plane-wave cutoff in Ry (default: 50)\n"
            "  ecutrho      Charge density cutoff in Ry (default: 8*ecutwfc)\n"
            "\n"
            "EXAMPLES\n"
            "  " << prog << " cif -pre si.cif 0.15\n"
            "  " << prog << " cif -pre al.cif 0.20 al.scf.in 40 320\n";
    }
    // ── dos ──────────────────────────────────────────────────────────────────
    else if (cmd == "dos") {
        std::cout <<
            "DESCRIPTION\n"
            "  Reads the DOS table produced by QE dos.x and generates PNG plots of the\n"
            "  density of states using Matplot++. The energy axis is shifted so that the\n"
            "  Fermi level sits at 0 eV. The Fermi energy can be supplied directly as a number\n"
            "  or extracted automatically from a QE SCF/NSCF output file.\n"
            "\n"
            "  If projwfc.x partial-DOS files (<prefix>.pdos_atm#*) are found in the same\n"
            "  directory, two additional plots are generated automatically:\n"
            "    <prefix>.pdos_elem.png  — elemental contributions (one curve per element)\n"
            "    <prefix>.pdos_orb.png   — orbital contributions summed (s, p, d, f)\n"
            "\n"
            "USAGE\n"
            "  " << prog << " dos -post <qe_dos_file> [fermi_eV|qe_output.out] [output_prefix]\n"
            "\n"
            "ARGUMENTS\n"
            "  qe_dos_file      Output of QE dos.x (e.g. si.dos)\n"
            "  fermi_eV         Fermi energy in eV, or path to QE scf/nscf output to extract it\n"
            "  output_prefix    Prefix for generated plot files (default: stem of dos file)\n"
            "\n"
            "EXAMPLES\n"
            "  " << prog << " dos -post si.dos 6.341\n"
            "  " << prog << " dos -post si.dos si.scf.out\n"
            "  " << prog << " dos -post si.dos si.scf.out si_plot\n"
            "\n"
            "NOTE\n"
            "  For PDOS, run 'projwfc.x' after 'dos.x' using the same prefix and outdir.\n"
            "  projwfc.x will create files like: <prefix>.pdos_atm#1(Si)_wfc#1(s)\n";
    }
    // ── band ─────────────────────────────────────────────────────────────────
    else if (cmd == "band") {
        if (sub == "-pre") {
            std::cout <<
                "DESCRIPTION\n"
                "  Automates the setup of a QE band-structure calculation. Given a CIF file\n"
                "  and an existing SCF input, it determines the crystal symmetry, selects the\n"
                "  standard high-symmetry k-path for that Bravais lattice, and writes two\n"
                "  ready-to-run input files: one for pw.x (non-SCF bands run) and one for\n"
                "  bands.x (post-processing). All settings (prefix, pseudo_dir, outdir,\n"
                "  cutoffs) are inherited from the SCF input.\n"
                "\n"
                "USAGE\n"
                "  " << prog << " band -pre <input.cif> <scf_input.in> [bands_pw.in] [bands_pp.in] [pts] [nbnd]\n"
                "\n"
                "ARGUMENTS\n"
                "  input.cif          CIF structure file (used to determine the k-path)\n"
                "  scf_input.in       Existing SCF input (prefix, pseudo_dir, outdir reused)\n"
                "  bands_pw.in        Output pw.x bands input (default: <stem>.bands.in)\n"
                "  bands_pp.in        Output bands.x input (default: <stem>.bands_pp.in)\n"
                "  pts                K-points per segment (default: 20)\n"
                "  nbnd               Total number of bands including empty bands\n"
                "                     (default: 0 = let QE decide; set e.g. 2x occupied to\n"
                "                     see conduction bands and the band gap)\n"
                "\n"
                "EXAMPLES\n"
                "  " << prog << " band -pre si.cif si.scf.in\n"
                "  " << prog << " band -pre si.cif si.scf.in si.bands.in si.bands_pp.in 40 8\n";
        } else if (sub == "-post") {
            std::cout <<
                "DESCRIPTION\n"
                "  Reads the band-structure data produced by QE bands.x (the .dat.gnu file)\n"
                "  and generates a PNG plot of the electronic band structure using Matplot++.\n"
                "  High-symmetry k-point positions are detected automatically from\n"
                "  discontinuities in the k-path; labels can optionally be supplied on the\n"
                "  command line. The energy axis is shifted to place the Fermi level at 0 eV.\n"
                "\n"
                "USAGE\n"
                "  " << prog << " band -post <qe_band_file> [fermi_eV|qe_output.out] [output_prefix] [labels]\n"
                "\n"
                "ARGUMENTS\n"
                "  qe_band_file     Output of QE bands.x (e.g. si.bands.dat.gnu)\n"
                "  fermi_eV         Fermi energy in eV, or path to QE scf/nscf output\n"
                "  output_prefix    Prefix for generated plot files (default: stem of band file)\n"
                "  labels           Comma-separated high-symmetry labels, e.g. L,G,X,W,K,G\n"
                "\n"
                "EXAMPLES\n"
                "  " << prog << " band -post si.bands.dat.gnu 6.341\n"
                "  " << prog << " band -post si.bands.dat.gnu si.scf.out si_bands L,G,X,W,K,G\n";
        } else if (sub == "-fat") {
            std::cout <<
                "DESCRIPTION\n"
                "  Plots a fat (projected) band structure. Thin gray lines show all bands;\n"
                "  colored scatter bubbles (area proportional to projection weight) show\n"
                "  the contribution of selected elements, atoms, or orbitals, computed from\n"
                "  the atomic_proj.xml produced by QE projwfc.x.\n"
                "\n"
                "  When band -pre is run it also writes a projwfc.x input file alongside\n"
                "  the bands_pp.in. After running pw.x (bands) and projwfc.x, call -fat.\n"
                "\n"
                "USAGE\n"
                "  " << prog << " band -fat <bands.dat.gnu> <atomic_proj.xml> [outprefix] [fermi_eV|qe_output.out] [filter ...]\n"
                "\n"
                "ARGUMENTS\n"
                "  bands.dat.gnu     Band-structure data file from bands.x (the .gnu file)\n"
                "  atomic_proj.xml   Projection data from projwfc.x (in <outdir>/<prefix>.save/)\n"
                "  outprefix         Output prefix for PNG files (default: stem of .gnu file)\n"
                "  fermi_eV          Fermi energy in eV, OR path to QE scf/nscf output\n"
                "  filter            One or more of (all combinable):\n"
                "    element=Si,Fe     filter by element name\n"
                "    atom=1,3          filter by 1-based atom index\n"
                "    orbital=s,p,d,f   filter by orbital type (s/p/d/f)\n"
                "\n"
                "GROUPING\n"
                "  element + orbital  → one PNG per (element, orbital) pair\n"
                "  element only       → one PNG per element\n"
                "  atom + orbital     → one PNG per (atom, orbital) pair\n"
                "  atom only          → one PNG per atom\n"
                "  orbital only       → one PNG per orbital, summed over all atoms\n"
                "  none               → auto: one PNG per distinct element\n"
                "  If >1 group, an additional combined PNG is produced.\n"
                "\n"
                "OUTPUT\n"
                "  <outprefix>.fatband_<group>.png  — one PNG per group\n"
                "  <outprefix>.fatband_combined.png — all groups overlaid (if >1 group)\n"
                "\n"
                "EXAMPLES\n"
                "  " << prog << " band -fat si.bands.dat.gnu tmp/si.save/atomic_proj.xml si 6.204\n"
                "  " << prog << " band -fat si.bands.dat.gnu tmp/si.save/atomic_proj.xml si 6.204 element=Si\n"
                "  " << prog << " band -fat si.bands.dat.gnu tmp/si.save/atomic_proj.xml si 6.204 orbital=s,p\n"
                "  " << prog << " band -fat si.bands.dat.gnu tmp/si.save/atomic_proj.xml si 6.204 element=Si orbital=s,p\n"
                "  " << prog << " band -fat bi2se3.bands.dat.gnu tmp/bi2se3.save/atomic_proj.xml bi2se3 12.5 element=Bi,Se\n"
                "  " << prog << " band -fat bi2se3.bands.dat.gnu tmp/bi2se3.save/atomic_proj.xml bi2se3 12.5 atom=1,2\n";
        } else {
            std::cout <<
                "USAGE\n"
                "  " << prog << " band -pre  <input.cif> <scf_input.in> [bands_pw.in] [bands_pp.in] [pts]\n"
                "  " << prog << " band -post <qe_band_file> [fermi_eV|qe_output.out] [output_prefix] [labels]\n"
                "  " << prog << " band -fat  <bands.dat.gnu> <atomic_proj.xml> [outprefix] [fermi_eV|out] [filter]\n"
                "\n"
                "Run '" << prog << " band -pre --help' or '" << prog << " band -post --help' or\n"
                "'" << prog << " band -fat --help' for details.\n";
        }
    }
    // ── kpath ─────────────────────────────────────────────────────────────────
    else if (cmd == "kpath") {
        std::cout <<
            "DESCRIPTION\n"
            "  Reads a CIF file, identifies the Bravais lattice and crystal family via spglib\n"
            "  (or ibrav fallback), and prints the standard high-symmetry k-path as a\n"
            "  K_POINTS crystal_b block ready to paste into a QE input file. Optionally writes\n"
            "  the block to a file.\n"
            "\n"
            "USAGE\n"
            "  " << prog << " kpath -pre <input.cif> [pts_per_segment] [output.kpath]\n"
            "\n"
            "ARGUMENTS\n"
            "  input.cif        CIF structure file\n"
            "  pts_per_segment  Points per k-path segment (default: 20)\n"
            "  output.kpath     Write K_POINTS block to this file (default: stdout only)\n"
            "\n"
            "EXAMPLES\n"
            "  " << prog << " kpath -pre si.cif\n"
            "  " << prog << " kpath -pre si.cif 40 si.kpath\n";
    }
    // ── elastic ───────────────────────────────────────────────────────────────
    else if (cmd == "elastic") {
        if (sub == "-pre") {
            std::cout <<
                "DESCRIPTION\n"
                "  Pre-processing step for computing elastic constants via the energy-strain\n"
                "  method. Reads an equilibrium QE SCF input, detects the crystal symmetry,\n"
                "  and generates a set of strained pw.x input files covering all\n"
                "  symmetry-independent strain patterns. Each pattern directory contains\n"
                "  inputs at evenly-spaced strain amplitudes from -max_delta to +max_delta.\n"
                "  After running all pw.x jobs, use 'elastic -post' to extract constants.\n"
                "\n"
                "USAGE\n"
                "  " << prog << " elastic -pre <scf_template.in> <outdir> [ndeltas] [max_delta]\n"
                "\n"
                "ARGUMENTS\n"
                "  scf_template.in  QE SCF input for the equilibrium structure\n"
                "                   (ibrav=0 with CELL_PARAMETERS, or ibrav+celldm)\n"
                "  outdir           Root directory where strain subdirectories are created\n"
                "  ndeltas          Number of strain points per pattern, odd (default: 7)\n"
                "  max_delta        Maximum strain magnitude, e.g. 0.04 = 4% (default: 0.04)\n"
                "\n"
                "OUTPUT STRUCTURE\n"
                "  outdir/\n"
                "    e1/   p0.0400/si.in  m0.0400/si.in  ...\n"
                "    e1e2/ p0.0400/si.in  ...\n"
                "    e4/   ...\n"
                "\n"
                "EXAMPLES\n"
                "  " << prog << " elastic -pre si.scf.in elastic_si\n"
                "  " << prog << " elastic -pre si.scf.in elastic_si 9 0.03\n"
                "\n"
                "NOTE\n"
                "  Crystal symmetry is detected automatically (spglib if available, else ibrav).\n"
                "  Only the symmetry-independent strain patterns are generated.\n";
        } else if (sub == "-post") {
            std::cout <<
                "DESCRIPTION\n"
                "  Post-processing step for elastic constants. Scans the output directories\n"
                "  created by 'elastic -pre', extracts the DFT total energies from each pw.x\n"
                "  output, fits a quadratic E(delta) polynomial per strain pattern, and\n"
                "  assembles the full 6x6 stiffness matrix (Cij, in GPa) using the appropriate\n"
                "  crystal-symmetry relations. Derived mechanical properties (bulk modulus K,\n"
                "  shear modulus G, Young's modulus E, Poisson's ratio nu) are reported via\n"
                "  the Voigt-Reuss-Hill average. A Born stability check is also performed.\n"
                "\n"
                "USAGE\n"
                "  " << prog << " elastic -post <scf_template.in> <outdir>\n"
                "\n"
                "ARGUMENTS\n"
                "  scf_template.in  Same template used with 'elastic -pre' (provides cell)\n"
                "  outdir           Same root directory used with 'elastic -pre'\n"
                "\n"
                "OUTPUT\n"
                "  Full 6x6 stiffness matrix (GPa), independent elastic constants, bulk\n"
                "  modulus, shear modulus, Young's modulus, Poisson's ratio, stability check.\n"
                "  Saves a report to outdir/elastic_results.txt\n"
                "\n"
                "EXAMPLES\n"
                "  " << prog << " elastic -post si.scf.in elastic_si\n";
        } else {
            std::cout <<
                "USAGE\n"
                "  " << prog << " elastic -pre  <scf_template.in> <outdir> [ndeltas] [max_delta]\n"
                "  " << prog << " elastic -post <scf_template.in> <outdir>\n"
                "\n"
                "WORKFLOW\n"
                "  1. " << prog << " elastic -pre  si.scf.in elastic_si\n"
                "  2. for d in elastic_si/*/*/; do (cd \"$d\" && pw.x < si.in > si.out); done\n"
                "  3. " << prog << " elastic -post si.scf.in elastic_si\n"
                "\n"
                "Run '" << prog << " elastic -pre --help' or '" << prog << " elastic -post --help' for details.\n";
        }
    }
    else if (cmd == "charge") {
        if (sub == "-pre") {
            std::cout <<
                "DESCRIPTION\n"
                "  Generates three pp.x input files for QE post-processing:\n"
                "    charge density         (plot_num=0)\n"
                "    charge-density difference (plot_num=6)\n"
                "    electron localisation function, ELF (plot_num=8)\n"
                "\n"
                "USAGE\n"
                "  " << prog << " charge -pre <scf.in> [outdir]\n"
                "\n"
                "ARGUMENTS\n"
                "  scf.in    QE SCF input file (used to extract prefix and outdir)\n"
                "  outdir    Directory where pp.x inputs are written (default: <stem>_charge/)\n"
                "\n"
                "OUTPUT\n"
                "  <outdir>/<prefix>.charge.pp.in\n"
                "  <outdir>/<prefix>.charge_diff.pp.in\n"
                "  <outdir>/<prefix>.elf.pp.in\n"
                "\n"
                "EXAMPLES\n"
                "  " << prog << " charge -pre si.scf.in si_charge\n";
        } else if (sub == "-post") {
            std::cout <<
                "DESCRIPTION\n"
                "  Reads a Gaussian cube file produced by QE pp.x and generates:\n"
                "    <outprefix>.slice.png  — three orthogonal midplane heatmaps (XY, XZ, YZ)\n"
                "    <outprefix>.3d.png     — 3-D scatter of high-value points (isosurface proxy)\n"
                "\n"
                "USAGE\n"
                "  " << prog << " charge -post <cube_file> [outprefix] [quantity]\n"
                "\n"
                "ARGUMENTS\n"
                "  cube_file   Cube file produced by pp.x (output_format=6)\n"
                "  outprefix   Output filename prefix (default: stem of cube_file)\n"
                "  quantity    One of: charge (default), charge_diff, elf\n"
                "              Affects colour map and thresholding\n"
                "\n"
                "EXAMPLES\n"
                "  " << prog << " charge -post si.charge.cube si_charge charge\n"
                "  " << prog << " charge -post si.charge_diff.cube si_charge_diff charge_diff\n"
                "  " << prog << " charge -post si.elf.cube si_elf elf\n";
        } else {
            std::cout <<
                "USAGE\n"
                "  " << prog << " charge -pre  <scf.in> [outdir]\n"
                "  " << prog << " charge -post <cube_file> [outprefix] [quantity]\n"
                "\n"
                "WORKFLOW\n"
                "  1. " << prog << " charge -pre si.scf.in si_charge\n"
                "  2. cd si_charge && pp.x < si.charge.pp.in > si.charge.pp.out\n"
                "  3. " << prog << " charge -post si.charge.cube si_charge charge\n"
                "\n"
                "Run '" << prog << " charge -pre --help' or '" << prog << " charge -post --help' for details.\n";
        }
    }
    else {
        std::cerr << "Unknown command '" << cmd << "'. Run '" << prog << " help' for the list.\n";
    }
}
