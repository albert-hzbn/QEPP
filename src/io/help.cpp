#include "qe/help.hpp"

#include <iostream>

void print_help(const char* prog) {
    std::cout <<
        "qepp — Quantum ESPRESSO post-processor and input generator\n"
        "\n"
        "USAGE\n"
        "  " << prog << " <command> -pre  [arguments]   (pre-processing)\n"
        "  " << prog << " <command> -run  [arguments]   (workflow execution)\n"
        "  " << prog << " <command> -post [arguments]   (post-processing)\n"
        "  " << prog << " help [command] [-pre|-run|-post]\n"
        "\n"
        "COMMANDS\n"
        "  cif     -pre   Generate a QE SCF input from a CIF file\n"
        "  dos     -post  Plot density of states from a QE dos.x output\n"
        "  band    -pre   Generate QE band-structure inputs (pw.x + bands.x) from a CIF\n"
        "  band    -post  Plot band structure from a QE bands.x output\n"
        "  kpath   -pre   Suggest a high-symmetry k-path from a CIF file\n"
        "  elastic -pre   Generate deformed SCF inputs for elastic constants\n"
        "  elastic -run   Run equilibrium + strained SCF jobs for elastic constants\n"
        "  elastic -post  Collect energies and compute elastic constants\n"
        "  charge  -pre   Generate pp.x inputs for charge density / charge diff / ELF\n"
        "  charge  -post  Plot volumetric data from a cube file\n"
        "  mag     -post  Tabulate per-atom magnetic moments from a QE pw.x output\n"
        "  stm     -pre   Generate pp.x input for STM simulation (ILDOS, Tersoff-Hamann)\n"
        "  stm     -post  Plot 2D constant-height STM map from a cube file\n"
        "  bader   -post  Tabulate per-atom Bader charges from ACF.dat\n"
        "  conv    -pre   Generate SCF inputs for ecutwfc or k-mesh convergence testing\n"
        "  conv    -post  Parse energies and plot convergence curve\n"
        "  struct  -post  Print lattice parameters, cell vectors, and atomic positions\n"
        "  parse   -post  Extract energy, Fermi level, forces, pressure, timing from QE output\n"
        "  qha     -pre   Generate scaled SCF inputs for quasi-harmonic approximation\n"
        "  qha     -post  Compute QHA thermal properties from phonopy yamls + DFT energies\n"
        "  qha_elastic -pre   Generate volume-scaled elastic + DFPT inputs\n"
        "  qha_elastic -run   Run the generated QHA elastic workflow\n"
        "  qha_elastic -post  Compute temperature-dependent elastic constants\n"
        "  phonon  -pre   Generate DFPT phonon inputs (ph.x, q2r.x, matdyn.x)\n"
        "  phonon  -post  Post-process DFPT matdyn.x output (DOS, band, HA)\n"
        "  phonon  -dos   Plot phonon DOS (phonopy or matdyn.x format)\n"
        "  phonon  -band  Plot phonon band structure (phonopy or matdyn.x format)\n"
        "  phonon  -ha    Harmonic approximation thermodynamics (phonopy or matdyn.x)\n"
        "\n"
        "Run '" << prog << " help <command>' or '" << prog << " <command> --help' for details.\n";
}

void print_help_command(const char* prog, const std::string& cmd,
                        const std::string& sub) {
    // ── cif ──────────────────────────────────────────────────────────────────
    if (cmd == "cif") {
        std::cout <<
            "USAGE\n"
            "  " << prog << " cif -pre <input.cif> <kspacing> [output.in] [ecutwfc] [ecutrho]\n"
            "\n"
            "ARGUMENTS\n"
            "  input.cif    CIF structure file\n"
            "  kspacing     K-point spacing in 1/Ang (e.g. 0.15); smaller = denser mesh\n"
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
            "USAGE\n"
            "  " << prog << " dos -post <qe_dos_file> [fermi_eV|qe_output.out] [output_prefix]\n"
            "\n"
            "ARGUMENTS\n"
            "  qe_dos_file      Output of QE dos.x (e.g. si.dos)\n"
            "  fermi_eV         Fermi energy in eV, or path to QE scf/nscf output\n"
            "  output_prefix    Output file prefix (default: stem of dos file)\n"
            "\n"
            "  If projwfc.x .pdos_atm#* files are present in the same directory,\n"
            "  elemental and orbital PDOS plots are generated automatically.\n"
            "\n"
            "EXAMPLES\n"
            "  " << prog << " dos -post si.dos 6.341\n"
            "  " << prog << " dos -post si.dos si.scf.out\n"
            "  " << prog << " dos -post si.dos si.scf.out si_plot\n";
    }
    // ── band ─────────────────────────────────────────────────────────────────
    else if (cmd == "band") {
        if (sub == "-pre") {
            std::cout <<
                "USAGE\n"
                "  " << prog << " band -pre <input.cif> <scf_input.in> [bands_pw.in] [bands_pp.in] [pts] [nbnd]\n"
                "\n"
                "ARGUMENTS\n"
                "  input.cif          CIF structure file (used to determine the k-path)\n"
                "  scf_input.in       Existing SCF input (prefix, pseudo_dir, outdir reused)\n"
                "  bands_pw.in        Output pw.x bands input (default: <stem>.bands.in)\n"
                "  bands_pp.in        Output bands.x input (default: <stem>.bands_pp.in)\n"
                "  pts                K-points per segment (default: 20)\n"
                "  nbnd               Number of bands (default: 0 = QE default)\n"
                "\n"
                "EXAMPLES\n"
                "  " << prog << " band -pre si.cif si.scf.in\n"
                "  " << prog << " band -pre si.cif si.scf.in si.bands.in si.bands_pp.in 40 8\n";
        } else if (sub == "-post") {
            std::cout <<
                "USAGE\n"
                "  " << prog << " band -post <qe_band_file> [fermi_eV|qe_output.out] [output_prefix] [labels]\n"
                "\n"
                "ARGUMENTS\n"
                "  qe_band_file     Output of QE bands.x (.dat.gnu file)\n"
                "  fermi_eV         Fermi energy in eV, or path to QE scf/nscf output\n"
                "  output_prefix    Output file prefix (default: stem of band file)\n"
                "  labels           Comma-separated high-symmetry labels, e.g. L,G,X,W,K,G\n"
                "\n"
                "EXAMPLES\n"
                "  " << prog << " band -post si.bands.dat.gnu 6.341\n"
                "  " << prog << " band -post si.bands.dat.gnu si.scf.out si_bands L,G,X,W,K,G\n";
        } else if (sub == "-fat") {
            std::cout <<
                "USAGE\n"
                "  " << prog << " band -fat <bands.dat.gnu> <atomic_proj.xml> [outprefix] [fermi_eV|qe_output.out] [filter ...]\n"
                "\n"
                "ARGUMENTS\n"
                "  bands.dat.gnu     Band-structure data file from bands.x\n"
                "  atomic_proj.xml   Projection data from projwfc.x (<outdir>/<prefix>.save/)\n"
                "  outprefix         Output prefix for PNG files (default: stem of .gnu file)\n"
                "  fermi_eV          Fermi energy in eV, or path to QE scf/nscf output\n"
                "  filter            element=Si,Fe | atom=1,3 | orbital=s,p,d,f  (combinable)\n"
                "\n"
                "EXAMPLES\n"
                "  " << prog << " band -fat si.bands.dat.gnu tmp/si.save/atomic_proj.xml si 6.204\n"
                "  " << prog << " band -fat si.bands.dat.gnu tmp/si.save/atomic_proj.xml si 6.204 element=Si orbital=s,p\n"
                "  " << prog << " band -fat bi2se3.bands.dat.gnu tmp/bi2se3.save/atomic_proj.xml bi2se3 12.5 element=Bi,Se\n";
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
                "USAGE\n"
                "  " << prog << " elastic -pre <scf_template.in> [outdir] [ndeltas] [max_delta]\n"
                "                     [--outdir D] [--ndeltas N] [--maxdelta D]\n"
                "\n"
                "ARGUMENTS\n"
                "  scf_template.in  Equilibrium QE SCF input (ibrav=0 or ibrav+celldm)\n"
                "  outdir           Output directory (default: <stem>_elastic)\n"
                "  ndeltas          Strain points per pattern, odd (default: 7)\n"
                "  max_delta        Max strain amplitude, e.g. 0.04 = 4% (default: 0.04)\n"
                "\n"
                "  Strain patterns are chosen automatically from crystal symmetry.\n"
                "\n"
                "EXAMPLES\n"
                "  " << prog << " elastic -pre si.scf.in\n"
                "  " << prog << " elastic -pre si.scf.in elastic_si 9 0.03\n"
                "  " << prog << " elastic -pre si.scf.in --outdir elastic_si --ndeltas 9 --maxdelta 0.03\n";
        } else if (sub == "-run") {
            std::cout <<
                "USAGE\n"
                "  " << prog << " elastic -run <scf_template.in> <outdir>\n"
                "  " << prog << " elastic -run <outdir_or_volume_dir>\n"
                "                     [--np N] [--ni N] [--nk N] [--nb N] [--nt N] [--nd N]\n"
                "\n"
                "  Runs the equilibrium SCF then all strained jobs. Skips stages with JOB DONE.\n"
                "  With a single directory argument, auto-detects the template and outdir.\n"
                "\n"
                "OPTIONS\n"
                "  --np N   MPI ranks (default: 1)\n"
                "  --ni N   QE -ni (image groups)\n"
                "  --nk N   QE -nk (k-point pools)\n"
                "  --nb N   QE -nb (band groups)\n"
                "  --nt N   QE -nt (task groups)\n"
                "  --nd N   QE -nd (diagonalization groups)\n"
                "\n"
                "EXAMPLES\n"
                "  " << prog << " elastic -run si.scf.in elastic_si --np 16 --nk 4\n"
                "  " << prog << " elastic -run tests/al_qha_el_litfix_new/v01 --np 16 --nk 4\n";
        } else if (sub == "-post") {
            std::cout <<
                "USAGE\n"
                "  " << prog << " elastic -post <scf_template.in> <outdir>\n"
                "  " << prog << " elastic -post <outdir_or_volume_dir>\n"
                "\n"
                "  Extracts DFT energies from each strained pw.x output, fits E(delta),\n"
                "  and assembles the full 6x6 Cij matrix (GPa) + VRH moduli.\n"
                "  Output saved to outdir/elastic_results.txt.\n"
                "\n"
                "EXAMPLES\n"
                "  " << prog << " elastic -post si.scf.in elastic_si\n"
                "  " << prog << " elastic -post tests/al_qha_el_litfix_new/v01\n";
        } else {
            std::cout <<
                "USAGE\n"
                "  " << prog << " elastic -pre  <scf_template.in> [outdir] [ndeltas] [max_delta] [--outdir D] [--ndeltas N] [--maxdelta D]\n"                "  " << prog << " elastic -run  <scf_template.in> <outdir> [--np N] [--ni N] [--nk N] [--nb N] [--nt N] [--nd N]\n"
                "  " << prog << " elastic -run  <outdir_or_volume_dir> [--np N] [--ni N] [--nk N] [--nb N] [--nt N] [--nd N]\n"
                "  " << prog << " elastic -post <scf_template.in> <outdir>\n"
                "  " << prog << " elastic -post <outdir_or_volume_dir>\n"
                "\n"
                "WORKFLOW\n"
                "  1. " << prog << " elastic -pre  si.scf.in elastic_si\n"
                "  2. " << prog << " elastic -run  si.scf.in elastic_si --np 16 --nk 4\n"
                "  3. " << prog << " elastic -post si.scf.in elastic_si\n"
                "\n"
                "Run '" << prog << " elastic -pre --help', '" << prog << " elastic -run --help', or '" << prog << " elastic -post --help' for details.\n";
        }
    }
    else if (cmd == "charge") {
        if (sub == "-pre") {
            std::cout <<
                "USAGE\n"
                "  " << prog << " charge -pre <scf.in> [outdir] [--outdir D]\n"
                "\n"
                "  Generates pp.x input files for charge density (plot_num=0),\n"
                "  charge-density difference (plot_num=6), and ELF (plot_num=8).\n"
                "  Output directory defaults to <stem>_charge.\n"
                "\n"
                "EXAMPLES\n"
                "  " << prog << " charge -pre si.scf.in\n"
                "  " << prog << " charge -pre si.scf.in si_charge\n";
        } else if (sub == "-post") {
            std::cout <<
                "USAGE\n"
                "  " << prog << " charge -post <cube_file> [outprefix] [quantity]\n"
                "\n"
                "ARGUMENTS\n"
                "  cube_file   Cube file produced by pp.x (output_format=6)\n"
                "  outprefix   Output filename prefix (default: stem of cube_file)\n"
                "  quantity    charge | charge_diff | elf (default: charge)\n"
                "\n"
                "EXAMPLES\n"
                "  " << prog << " charge -post si.charge.cube\n"
                "  " << prog << " charge -post si.charge_diff.cube si_charge charge_diff\n";
        } else {
            std::cout <<
                "USAGE\n"
                "  " << prog << " charge -pre  <scf.in> [outdir] [--outdir D]\n"
                "  " << prog << " charge -post <cube_file> [outprefix] [quantity]\n"
                "\n"
                "WORKFLOW\n"
                "  1. " << prog << " charge -pre si.scf.in si_charge\n"
                "  2. cd si_charge && pp.x < si.charge.pp.in > si.charge.pp.out\n"
                "  3. " << prog << " charge -post si.charge.cube\n"
                "\n"
                "Run '" << prog << " charge -pre --help' or '" << prog << " charge -post --help' for details.\n";
        }
    }
    else if (cmd == "mag") {
        std::cout <<
            "USAGE\n"
            "  " << prog << " mag -post <qe_output.out> [output_prefix]\n"
            "\n"
            "  Prints per-atom magnetic moments from a QE pw.x output.\n"
            "  Supports collinear (nspin=2) and non-collinear calculations.\n"
            "\n"
            "EXAMPLES\n"
            "  " << prog << " mag -post fe.scf.out\n"
            "  " << prog << " mag -post fe3o4.scf.out fe3o4_mag\n";
    }
    else if (cmd == "stm") {
        if (sub == "-pre") {
            std::cout <<
                "USAGE\n"
                "  " << prog << " stm -pre <scf.in> [bias_eV] [outdir] [--outdir D]\n"
                "\n"
                "  Generates a pp.x input for STM simulation (ILDOS, plot_num=5).\n"
                "  Positive bias probes empty states [0, bias]; negative probes filled [bias, 0].\n"
                "  Output directory defaults to <stem>_stm.\n"
                "\n"
                "EXAMPLES\n"
                "  " << prog << " stm -pre si.scf.in\n"
                "  " << prog << " stm -pre si.scf.in 0.5 si_stm\n"
                "  " << prog << " stm -pre si.scf.in --outdir si_stm\n";
        } else if (sub == "-post") {
            std::cout <<
                "USAGE\n"
                "  " << prog << " stm -post <stm.cube> [output_prefix] [height_ang]\n"
                "\n"
                "ARGUMENTS\n"
                "  stm.cube       Cube file from pp.x STM run\n"
                "  output_prefix  Output prefix for PNG (default: stem of cube file)\n"
                "  height_ang     Tip height in Ang above cube origin (-1 = midplane, default)\n"
                "\n"
                "EXAMPLES\n"
                "  " << prog << " stm -post si.stm.cube\n"
                "  " << prog << " stm -post si.stm.cube si_stm 12.5\n";
        } else {
            std::cout <<
                "USAGE\n"
                "  " << prog << " stm -pre  <scf.in> [bias_eV] [outdir] [--outdir D]\n"
                "  " << prog << " stm -post <stm.cube> [output_prefix] [height_ang]\n"
                "\n"
                "WORKFLOW\n"
                "  1. " << prog << " stm -pre si.scf.in 1.0 si_stm\n"
                "  2. cd si_stm && pp.x < si.stm.pp.in > si.stm.pp.out\n"
                "  3. " << prog << " stm -post si.stm.cube\n"
                "\n"
                "Run '" << prog << " stm -pre --help' or '" << prog << " stm -post --help' for details.\n";
        }
    }
    else if (cmd == "bader") {
        std::cout <<
            "USAGE\n"
            "  " << prog << " bader -post <ACF.dat> [scf.in] [output_prefix]\n"
            "\n"
            "ARGUMENTS\n"
            "  ACF.dat        Atomic charge file from Henkelman bader\n"
            "  scf.in         QE SCF input for element labels (optional but recommended)\n"
            "  output_prefix  Output prefix (default: stem of ACF.dat)\n"
            "\n"
            "WORKFLOW\n"
            "  1. qepp charge -pre si.scf.in si_charge\n"
            "  2. pp.x < si_charge/si.charge.pp.in > si_charge/si.charge.pp.out\n"
            "  3. bader si.charge.cube   (Henkelman bader binary)\n"
            "  4. qepp bader -post ACF.dat si.scf.in\n"
            "\n"
            "EXAMPLES\n"
            "  " << prog << " bader -post ACF.dat\n"
            "  " << prog << " bader -post ACF.dat si.scf.in si_bader\n";
    }
    // ── conv ─────────────────────────────────────────────────────────────────
    else if (cmd == "conv") {
        if (sub == "-pre") {
            std::cout <<
                "USAGE\n"
                "  " << prog << " conv -pre <scf.in> <ecutwfc|kspacing> <min> <max> <step> [outdir] [--outdir D]\n"
                "\n"
                "ARGUMENTS\n"
                "  scf.in     Template QE SCF input\n"
                "  ecutwfc    Vary plane-wave cutoff (Ry)\n"
                "  kspacing   Vary k-point spacing (1/Ang); requires CELL_PARAMETERS\n"
                "  min/max/step  Parameter range (inclusive)\n"
                "  outdir     Output directory (default: conv_ecutwfc or conv_kspacing)\n"
                "\n"
                "EXAMPLES\n"
                "  " << prog << " conv -pre si.scf.in ecutwfc 20 80 10\n"
                "  " << prog << " conv -pre si.scf.in kspacing 0.05 0.30 0.05 si_kconv\n"
                "  " << prog << " conv -pre si.scf.in ecutwfc 20 80 10 --outdir si_ecut\n";
        } else if (sub == "-post") {
            std::cout <<
                "USAGE\n"
                "  " << prog << " conv -post <outdir> <ecutwfc|kspacing> [output_prefix]\n"
                "\n"
                "ARGUMENTS\n"
                "  outdir            Directory created by conv -pre\n"
                "  ecutwfc|kspacing  Must match what was used in -pre\n"
                "  output_prefix     Output prefix (default: outdir/param)\n"
                "\n"
                "EXAMPLES\n"
                "  " << prog << " conv -post conv_ecutwfc ecutwfc\n"
                "  " << prog << " conv -post si_kconv kspacing\n";
        } else {
            std::cout <<
                "USAGE\n"
                "  " << prog << " conv -pre  <scf.in> <ecutwfc|kspacing> <min> <max> <step> [outdir] [--outdir D]\n"
                "  " << prog << " conv -post <outdir> <ecutwfc|kspacing> [output_prefix]\n"
                "\n"
                "WORKFLOW\n"
                "  1. " << prog << " conv -pre si.scf.in ecutwfc 20 80 10 si_ecut\n"
                "  2. for d in si_ecut/*/; do (cd \"$d\" && pw.x < scf.in > scf.out); done\n"
                "  3. " << prog << " conv -post si_ecut ecutwfc\n"
                "\n"
                "Run '" << prog << " conv -pre --help' or '" << prog << " conv -post --help' for details.\n";
        }
    }
    // ── struct ────────────────────────────────────────────────────────────────
    else if (cmd == "struct") {
        std::cout <<
            "USAGE\n"
            "  " << prog << " struct -post <structure_file> [output_prefix]\n"
            "                       [--sro] [--nshells N] [--tol T] [--source input|output|auto]\n"
            "\n"
            "ARGUMENTS\n"
            "  structure_file  QE pw.x SCF input or QE output with final coordinates\n"
            "  output_prefix   Output prefix (default: stem of structure_file)\n"
            "  --sro           Compute Warren-Cowley short-range-order parameters\n"
            "  --nshells N     Number of neighbor shells for SRO (default: 2)\n"
            "  --tol T         Shell clustering tolerance in Ang (default: 1e-3)\n"
            "  --source        input | output | auto (default: auto)\n"
            "\n"
            "EXAMPLES\n"
            "  " << prog << " struct -post si.scf.in\n"
            "  " << prog << " struct -post ni_co.scf.in ni_co --sro --nshells 3\n"
            "  " << prog << " struct -post ni_co.relax.out --sro --source output\n";
    }
    // ── parse ─────────────────────────────────────────────────────────────────
    else if (cmd == "parse") {
        std::cout <<
            "USAGE\n"
            "  " << prog << " parse -post <qe.out> [output_prefix]\n"
            "\n"
            "  Extracts total energy, Fermi energy, forces, pressure, SCF status,\n"
            "  and wall time from a QE pw.x output. Saves <prefix>.parse.txt.\n"
            "\n"
            "EXAMPLES\n"
            "  " << prog << " parse -post si.scf.out\n"
            "  " << prog << " parse -post al.relax.out al_relax\n";
    }
    // ── qha ──────────────────────────────────────────────────────────────────
    else if (cmd == "qha") {
        if (sub.empty() || sub == "-pre") {
            std::cout <<
                "USAGE\n"
                "  " << prog << " qha -pre <scf.in> [outdir] [--nvolumes N] [--range R] [--outdir D]\n"
                "\n"
                "ARGUMENTS\n"
                "  scf.in        Equilibrium QE SCF input (must have CELL_PARAMETERS angstrom)\n"
                "  outdir        Output directory (default: <stem>_qha)\n"
                "  --nvolumes N  Number of volume points (default: 7, minimum: 4)\n"
                "  --range R     Volume range in % (default: 10 → ±5%)\n"
                "  --outdir D    Same as positional outdir (flag form)\n"
                "\n"
                "EXAMPLES\n"
                "  " << prog << " qha -pre si.scf.in\n"
                "  " << prog << " qha -pre si.scf.in si_qha --nvolumes 9 --range 12\n"
                "  " << prog << " qha -pre si.scf.in --nvolumes 9 --range 12 --outdir si_qha\n";
        }
        if (sub.empty() || sub == "-post") {
            std::cout <<
                "USAGE\n"
                "  " << prog << " qha -post <qha_summary.in> [output_prefix]\n"
                "\n"
                "  Reads a summary file (volume, energy, yaml_path per line)\n"
                "  and computes QHA thermal properties (V(T), G(T), B_T(T), alpha(T), Cv, Cp).\n"
                "\n"
                "EXAMPLES\n"
                "  " << prog << " qha -post si_qha/qha_summary.in\n"
                "  " << prog << " qha -post si_qha/qha_summary.in si_qha_results\n";
        }
    }
    // ── phonon ───────────────────────────────────────────────────────────────
    else if (cmd == "phonon") {
        if (sub.empty() || sub == "-pre") {
            std::cout <<
                "USAGE\n"
                "  " << prog << " phonon -pre <scf.in> [outdir] [--outdir D]\n"
                "                         [--nq NQ1 NQ2 NQ3] [--nq-dos NQ1 NQ2 NQ3]\n"
                "                         [--epsil] [--asr TYPE]\n"
                "\n"
                "ARGUMENTS\n"
                "  scf.in            pw.x SCF input\n"
                "  outdir            Output directory (default: <stem>_phonon)\n"
                "  --nq NQ1 NQ2 NQ3  q-mesh for ph.x (default: 4 4 4)\n"
                "  --nq-dos N N N    q-mesh for matdyn DOS (default: 8 8 8)\n"
                "  --epsil           Include Born charges / dielectric constant\n"
                "  --asr TYPE        Acoustic sum rule: no/simple/crystal (default: simple)\n"
                "\n"
                "EXAMPLES\n"
                "  " << prog << " phonon -pre si_scf.in\n"
                "  " << prog << " phonon -pre si_scf.in --nq 6 6 6 ph_inputs/\n"
                "  " << prog << " phonon -pre si_scf.in --outdir ph_inputs/ --nq 6 6 6\n";
        }
        if (sub.empty() || sub == "-post") {
            std::cout <<
                "USAGE\n"
                "  " << prog << " phonon -post <prefix> [--dos] [--band] [--ha]\n"
                "                          [--natom N] [--tmin T] [--tmax T] [--dt T]\n"
                "                          [--labels G,X,M,G] [--qlabels <file>]\n"
                "                          [output_prefix]\n"
                "\n"
                "ARGUMENTS\n"
                "  prefix            Stem used in matdyn_dos.in / matdyn_band.in\n"
                "  --dos/--band/--ha  Select outputs (default: all three)\n"
                "  --natom N         Atoms per cell for per-atom normalisation (default: 1)\n"
                "  --tmin/tmax/dt    Temperature range in K (defaults: 50/1000/50)\n"
                "  --labels G,X,...  High-symmetry point labels (comma-separated)\n"
                "  --qlabels file    File with 'label kx ky kz' lines\n"
                "  output_prefix     Output prefix (default: same as prefix)\n"
                "\n"
                "EXAMPLES\n"
                "  " << prog << " phonon -post si\n"
                "  " << prog << " phonon -post si --ha --natom 2 --tmax 1200\n"
                "  " << prog << " phonon -post si --band --labels G,X,W,L,G si_bands\n";
        }
        if (sub.empty() || sub == "-dos") {
            std::cout <<
                "USAGE\n"
                "  " << prog << " phonon -dos <dos_file> [--pdos projected_dos.dat]\n"
                "                             [--labels a,b,c] [output_prefix]\n"
                "\n"
                "ARGUMENTS\n"
                "  dos_file           phonopy total_dos.dat or matdyn.x .dos file\n"
                "  --pdos file        phonopy projected_dos.dat (phonopy workflow only)\n"
                "  --labels a,b,...   Atom labels for PDOS legend\n"
                "  output_prefix      Output prefix (default: stem of input)\n"
                "\n"
                "EXAMPLES\n"
                "  " << prog << " phonon -dos total_dos.dat\n"
                "  " << prog << " phonon -dos si.phonon.dos si_dos\n";
        }
        if (sub.empty() || sub == "-band") {
            std::cout <<
                "USAGE\n"
                "  " << prog << " phonon -band <file> [--labels G,X,M,G]\n"
                "                               [--qlabels <file>] [output_prefix]\n"
                "\n"
                "ARGUMENTS\n"
                "  file               phonopy band.yaml or matdyn.x .freq file\n"
                "  --labels G,X,...   High-symmetry point labels (comma-separated)\n"
                "  --qlabels file     File with 'label kx ky kz' lines (matdyn workflow)\n"
                "  output_prefix      Output prefix (default: stem of input)\n"
                "\n"
                "EXAMPLES\n"
                "  " << prog << " phonon -band band.yaml\n"
                "  " << prog << " phonon -band si.phonon_band.freq --labels G,X,W,L,G\n";
        }
        if (sub.empty() || sub == "-ha") {
            std::cout <<
                "USAGE\n"
                "  " << prog << " phonon -ha <dos_file> [--tmin T] [--tmax T] [--dt T]\n"
                "                                       [--natom N] [output_prefix]\n"
                "\n"
                "ARGUMENTS\n"
                "  dos_file       phonopy total_dos.dat or matdyn.x .dos file\n"
                "  --tmin/tmax/dt Temperature range in K (defaults: 50/1000/50)\n"
                "  --natom N      Atoms per cell for per-atom normalisation (default: 1)\n"
                "  output_prefix  Output prefix (default: stem of input)\n"
                "\n"
                "EXAMPLES\n"
                "  " << prog << " phonon -ha total_dos.dat --natom 2\n"
                "  " << prog << " phonon -ha si.phonon.dos --natom 2 --tmax 1500 si_ha\n";
        }
    }
    // ── qha_elastic ──────────────────────────────────────────────────────────
    else if (cmd == "qha_elastic") {
        if (sub.empty() || sub == "-pre") {
            std::cout <<
                "USAGE\n"
                "  " << prog << " qha_elastic -pre <scf.in> [outdir]\n"
                "                    [--nvolumes N] [--range R]\n"
                "                    [--ndeltas N]  [--maxdelta D]\n"
                "                    [--nq N]       [--nq-dos N]\n"
                "                    [--tr2ph X]    [--outdir D]\n"
                "\n"
                "OPTIONS\n"
                "  --nvolumes N    Number of volumes (default: 7)\n"
                "  --range R       Volume range in % (default: 10 → ±5%)\n"
                "  --ndeltas N     Strain points per pattern, odd (default: 7)\n"
                "  --maxdelta D    Max strain amplitude (default: 0.04)\n"
                "  --nq N          Isotropic DFPT q-mesh (default: 4)\n"
                "  --nq-dos N      DOS mesh for matdyn.x (default: 16)\n"
                "  --tr2ph X       ph.x convergence threshold (default: 1e-14)\n"
                "  --outdir D      Output directory (default: <stem>_qha_elastic)\n"
                "\n"
                "EXAMPLES\n"
                "  " << prog << " qha_elastic -pre si_scf.in\n"
                "  " << prog << " qha_elastic -pre si_scf.in si_qha_el\n"
                "  " << prog << " qha_elastic -pre si_scf.in --nvolumes 9 --nq 2 --outdir si_qha_el\n";
        }
        if (sub.empty() || sub == "-run") {
            std::cout <<
                "USAGE\n"
                "  " << prog << " qha_elastic -run [dataset_dir]\n"
                "                     [--np N] [--ni N] [--nk N] [--nb N] [--nt N] [--nd N]\n"
                "                     [--exclude v04,v07]\n"
                "\n"
                "  Runs all volumes in the dataset (SCF, elastic SCFs, ph.x, q2r.x, matdyn.x).\n"
                "  Skips stages already containing JOB DONE. Defaults to current directory.\n"
                "\n"
                "OPTIONS\n"
                "  --np N          MPI ranks (default: 1)\n"
                "  --ni/nk/nb/nt/nd N  QE parallelization flags\n"
                "  --exclude LIST  Comma-separated volume dirs to skip\n"
                "\n"
                "EXAMPLES\n"
                "  " << prog << " qha_elastic -run si_qha_el --np 20 --nk 4\n"
                "  " << prog << " qha_elastic -run --np 20 --nk 4 --exclude v04\n";
        }
        if (sub.empty() || sub == "-post") {
            std::cout <<
                "USAGE\n"
                "  " << prog << " qha_elastic -post [qha_elastic_summary.in|dataset_dir]\n"
                "                     [output_prefix]\n"
                "                     [--tmin T] [--tmax T] [--dt T] [--exclude v04]\n"
                "\n"
                "  Reads elastic results and phonon DOS per volume, then computes\n"
                "  temperature-dependent C_ij(T) tensors and VRH moduli vs temperature.\n"
                "  Static energies are read from qe.out if the energy column is not numeric.\n"
                "\n"
                "OPTIONS\n"
                "  --tmin/tmax/dt  Temperature range in K (defaults: 0/1500/10)\n"
                "  --exclude       Comma-separated volume dirs to ignore\n"
                "\n"
                "EXAMPLES\n"
                "  " << prog << " qha_elastic -post si_qha_el\n"
                "  " << prog << " qha_elastic -post --exclude v04 --tmax 1000\n";
        }
    }
    else {
        std::cerr << "Unknown command '" << cmd << "'. Run '" << prog << " help' for the list.\n";
    }
}
