#include <Eigen/Dense>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "qe/band.hpp"
#include "qe/cif.hpp"
#include "qe/dos.hpp"
#include "qe/elastic.hpp"
#include "qe/qe_input.hpp"
#include "qe/types.hpp"
#include "qe/utils.hpp"

// Forward declarations so handlers inside namespace qe can call these.
static void print_help(const char* prog);
static void print_help_command(const char* prog, const std::string& cmd,
                               const std::string& sub = "");

namespace qe {

// s = argv-start offset: 0 when called as "qepp cif ...", 1 when called as "qepp cif -pre ..."
static int handle_cif_mode(int argc, char** argv, int s = 0) {
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
    const std::string outPath  = (argc >= 5 + s) ? argv[4 + s] : defaultOutput;
    const int ecutwfc          = (argc >= 6 + s) ? std::stoi(argv[5 + s]) : 50;
    const int ecutrho          = (argc >= 7 + s) ? std::stoi(argv[6 + s]) : (8 * ecutwfc);

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

static int handle_dos_mode(int argc, char** argv, int s = 0) {
    if (argc < 3 + s || argc > 5 + s) {
        print_help_command(argv[0], "dos", "-post");
        return 1;
    }

    const std::string dosPath = argv[2 + s];
    double fermiEv = 0.0;
    std::string outPrefix = stem_from_path(dosPath);

    if (argc >= 4 + s) {
        const std::string arg = argv[3 + s];
        if (!try_parse_double(arg, fermiEv)) {
            fermiEv = extract_fermi_from_qe_output(arg);
            std::cout << "Extracted Fermi energy from QE output: " << std::fixed
                      << std::setprecision(6) << fermiEv << " eV\n";
        }
    }
    if (argc >= 5 + s) {
        outPrefix = argv[4 + s];
    }

    const auto dosRows = parse_dos_table(dosPath);
    write_dos_plot_bundle(dosPath, dosRows, fermiEv, outPrefix);
    write_pdos_plots(dosPath, dosRows, fermiEv, outPrefix);
    return 0;
}

static int handle_band_post_mode(int argc, char** argv, int s = 0) {
    if (argc < 3 + s || argc > 6 + s) {
        print_help_command(argv[0], "band", "-post");
        return 1;
    }

    const std::string bandPath = argv[2 + s];
    double fermiEv = 0.0;
    std::string outPrefix = stem_from_path(bandPath);

    if (argc >= 4 + s) {
        const std::string arg = argv[3 + s];
        if (!try_parse_double(arg, fermiEv)) {
            fermiEv = extract_fermi_from_qe_output(arg);
            std::cout << "Extracted Fermi energy from QE output: " << std::fixed
                      << std::setprecision(6) << fermiEv << " eV\n";
        }
    }
    if (argc >= 5 + s) {
        outPrefix = argv[4 + s];
    }

    auto bandData = parse_band_table(bandPath);

    if (argc >= 6 + s) {
        const std::string labelsArg = argv[5 + s];
        std::vector<std::string> labels;
        std::istringstream ss(labelsArg);
        std::string token;
        while (std::getline(ss, token, ',')) {
            labels.push_back(trim(token));
        }
        if (labels.size() == bandData.kLabelMarks.size()) {
            for (size_t i = 0; i < labels.size(); ++i) {
                bandData.kLabelMarks[i].second = labels[i];
            }
        } else {
            std::cerr << "Warning: " << labels.size() << " labels given but "
                      << bandData.kLabelMarks.size()
                      << " high-symmetry points detected; labels ignored.\n";
        }
    }

    write_band_plot_bundle(bandPath, bandData, fermiEv, outPrefix);
    return 0;
}

static int handle_kpath_mode(int argc, char** argv, int s = 0) {
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

static int handle_band_pre_mode(int argc, char** argv, int s = 0) {
    if (argc < 4 + s || argc > 7 + s) {
        print_help_command(argv[0], "band", "-pre");
        return 1;
    }

    const std::string cifPath      = argv[2 + s];
    const std::string scfInputPath = argv[3 + s];
    const std::string defaultPrefix = stem_from_path(cifPath);
    const std::string bandsPwPath  = (argc >= 5 + s) ? argv[4 + s] : (defaultPrefix + ".bands.in");
    const std::string bandsPpPath  = (argc >= 6 + s) ? argv[5 + s] : (defaultPrefix + ".bands_pp.in");
    const int pointsPerSegment     = (argc >= 7 + s) ? std::stoi(argv[6 + s]) : 20;

    const CifStructure structure = parse_cif(cifPath);
    const SymmetryKPath kpath = suggest_kpath_from_cif(structure);

    write_bands_input_from_scf_template(scfInputPath, bandsPwPath, kpath,
                                        pointsPerSegment);
    const std::string filbandName = defaultPrefix + ".bands.dat";
    write_bands_pp_input_from_scf_template(scfInputPath, bandsPpPath, filbandName);

    std::cout << "Generated QE bands pw.x input: " << bandsPwPath << "\n";
    std::cout << "Generated QE bands.x input: " << bandsPpPath << "\n";
    std::cout << "Suggested k-path family: " << kpath.family << "\n";
    std::cout << "Path:";
    for (size_t i = 0; i < kpath.nodes.size(); ++i) {
        std::cout << " " << kpath.nodes[i].label;
        if (i + 1 < kpath.nodes.size()) std::cout << "-";
    }
    std::cout << "\n";
    std::cout << "Expected filband output name: " << filbandName << "\n";
    return 0;
}

static int handle_elastic_pre_mode(int argc, char** argv, int s = 0) {
    if (argc < 4 + s) {
        print_help_command(argv[0], "elastic", "-pre");
        return 1;
    }
    const std::string tmpl   = argv[2 + s];
    const std::string outDir = argv[3 + s];
    const int nDeltas        = (argc >= 5 + s) ? std::stoi(argv[4 + s]) : 7;
    const double maxDelta    = (argc >= 6 + s) ? std::stod(argv[5 + s]) : 0.04;
    generate_elastic_inputs(tmpl, outDir, nDeltas, maxDelta);
    return 0;
}

static int handle_elastic_post_mode(int argc, char** argv, int s = 0) {
    if (argc < 4 + s) {
        print_help_command(argv[0], "elastic", "-post");
        return 1;
    }
    const std::string tmpl   = argv[2 + s];
    const std::string outDir = argv[3 + s];
    const ElasticResults res = compute_elastic_properties(tmpl, outDir);
    print_elastic_report(res, res.crystalFamily, outDir + "/elastic_results.txt");
    return 0;
}

}  // namespace qe

static void print_help(const char* prog) {
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
        "\n"
        "Run '" << prog << " help <command>' or '" << prog << " <command> --help' for details.\n";
}

static void print_help_command(const char* prog, const std::string& cmd,
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
                "  " << prog << " band -pre <input.cif> <scf_input.in> [bands_pw.in] [bands_pp.in] [pts]\n"
                "\n"
                "ARGUMENTS\n"
                "  input.cif          CIF structure file (used to determine the k-path)\n"
                "  scf_input.in       Existing SCF input (prefix, pseudo_dir, outdir reused)\n"
                "  bands_pw.in        Output pw.x bands input (default: <stem>.bands.in)\n"
                "  bands_pp.in        Output bands.x input (default: <stem>.bands_pp.in)\n"
                "  pts                K-points per segment (default: 20)\n"
                "\n"
                "EXAMPLES\n"
                "  " << prog << " band -pre si.cif si.scf.in\n"
                "  " << prog << " band -pre si.cif si.scf.in si.bands.in si.bands_pp.in 40\n";
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
        } else {
            std::cout <<
                "USAGE\n"
                "  " << prog << " band -pre  <input.cif> <scf_input.in> [bands_pw.in] [bands_pp.in] [pts]\n"
                "  " << prog << " band -post <qe_band_file> [fermi_eV|qe_output.out] [output_prefix] [labels]\n"
                "\n"
                "Run '" << prog << " band -pre --help' or '" << prog << " band -post --help' for details.\n";
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
    else {
        std::cerr << "Unknown command '" << cmd << "'. Run '" << prog << " help' for the list.\n";
    }
}

int main(int argc, char** argv) {
    try {
        if (argc < 2) {
            print_help(argv[0]);
            return 1;
        }

        const std::string mode = qe::to_lower(argv[1]);

        // qepp help [cmd] [-pre|-post]  /  qepp --help [cmd]  /  qepp -h [cmd]
        if (mode == "help" || mode == "--help" || mode == "-h") {
            const std::string cmd = (argc >= 3) ? qe::to_lower(argv[2]) : "";
            const std::string sub = (argc >= 4) ? qe::to_lower(argv[3]) : "";
            if (cmd.empty()) { print_help(argv[0]); return 0; }
            print_help_command(argv[0], cmd, sub);
            return 0;
        }

        // qepp <cmd> --help  /  qepp <cmd> -pre --help  /  qepp <cmd> -post --help
        if (argc >= 3) {
            const std::string a2 = qe::to_lower(argv[2]);
            if (a2 == "--help" || a2 == "-h" || a2 == "help") {
                print_help_command(argv[0], mode, "");
                return 0;
            }
            if (argc >= 4) {
                const std::string a3 = qe::to_lower(argv[3]);
                if (a3 == "--help" || a3 == "-h" || a3 == "help") {
                    print_help_command(argv[0], mode, a2);
                    return 0;
                }
            }
        }

        // Resolve sub-flag (-pre / -post)
        const std::string sub = (argc >= 3) ? qe::to_lower(argv[2]) : "";

        if (mode == "cif") {
            if (sub != "-pre") { print_help_command(argv[0], "cif", ""); return 1; }
            return qe::handle_cif_mode(argc, argv, 1);
        }
        if (mode == "dos") {
            if (sub != "-post") { print_help_command(argv[0], "dos", ""); return 1; }
            return qe::handle_dos_mode(argc, argv, 1);
        }
        if (mode == "band") {
            if (sub == "-pre")  return qe::handle_band_pre_mode(argc, argv, 1);
            if (sub == "-post") return qe::handle_band_post_mode(argc, argv, 1);
            print_help_command(argv[0], "band", "");
            return 1;
        }
        if (mode == "kpath") {
            if (sub != "-pre") { print_help_command(argv[0], "kpath", ""); return 1; }
            return qe::handle_kpath_mode(argc, argv, 1);
        }
        if (mode == "elastic") {
            if (sub == "-pre")  return qe::handle_elastic_pre_mode(argc, argv, 1);
            if (sub == "-post") return qe::handle_elastic_post_mode(argc, argv, 1);
            print_help_command(argv[0], "elastic", "");
            return 1;
        }

        throw std::runtime_error("Unknown command: '" + mode +
                                 "'. Run '" + std::string(argv[0]) + " help' for usage.");
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << "\n";
        return 1;
    }
}


