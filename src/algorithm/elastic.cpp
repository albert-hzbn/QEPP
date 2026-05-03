#include "qe/elastic.hpp"
#include "qe/utils.hpp"

#ifdef HAVE_SPGLIB
#  include <spglib.h>
#endif

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace fs = std::filesystem;

namespace qe {

// ─────────────────────────────────────────────────────────────────────────────
//  Physical constants
// ─────────────────────────────────────────────────────────────────────────────
static constexpr double kB       = 1.380649e-23;    // J/K
static constexpr double hbar     = 1.054571817e-34; // J·s
static constexpr double GPa2Pa   = 1.0e9;
static constexpr double Ry2eV    = 13.605693122;
static constexpr double eV2J     = 1.602176634e-19;
static constexpr double Ang2m    = 1.0e-10;
static constexpr double bohr2ang = 0.529177210903;
static constexpr double amuToKg  = 1.66053906660e-27;

double extract_total_energy(const std::string& path);

// ─────────────────────────────────────────────────────────────────────────────
//  Crystal family from space-group number (ITA 2016)
// ─────────────────────────────────────────────────────────────────────────────
static std::string crystal_family_from_sg(int sg) {
    if (sg < 1 || sg > 230) return "triclinic";
    if (sg <= 2)   return "triclinic";
    if (sg <= 15)  return "monoclinic";
    if (sg <= 74)  return "orthorhombic";
    if (sg <= 88)  return "tetragonal_I";   // C4, S4, C4h: 7 constants (incl. C16)
    if (sg <= 142) return "tetragonal";     // D4, C4v, D2d, D4h: 6 constants
    if (sg <= 148) return "trigonal_I";     // C3, S6: 7 constants (incl. C15)
    if (sg <= 167) return "trigonal";       // D3, C3v, D3d: 6 constants
    if (sg <= 194) return "hexagonal";      // 5 constants
    return "cubic";                         // 3 constants
}

// ─────────────────────────────────────────────────────────────────────────────
//  QEInput struct
// ─────────────────────────────────────────────────────────────────────────────
struct QEInput {
    Eigen::Matrix3d cell = Eigen::Matrix3d::Zero();
    std::vector<std::pair<std::string, Eigen::Vector3d>> atoms;
    std::vector<std::string> speciesLines;
    std::string prefix    = "pwscf";
    std::string pseudoDir = "./pseudo";
    std::string outdir    = "./tmp";
    int ecutwfc = 50, ecutrho = 400;
    std::string kpointsBlock;
    double volume = 0.0;
    int ibrav = 0;
    std::array<double,6> celldm = {0,0,0,0,0,0};
    // smearing / occupations (metals)
    std::string occupations;  // e.g. "smearing"
    std::string smearing;     // e.g. "mp"
    double degauss = 0.0;
};

static QEInput parse_qe_input(const std::string& path) {
    const auto lines = load_lines(path);
    QEInput qi;
    bool inCell = false, inAtoms = false, inSpecies = false, inKpoints = false;
    int cellRow = 0;
    std::ostringstream kbuf;

    auto stripVal = [](const std::string& t, size_t eq) -> std::string {
        std::string v = trim(t.substr(eq + 1));
        while (!v.empty() && (v.back() == ',' || v.back() == ' ')) v.pop_back();
        return v;
    };

    for (const auto& raw : lines) {
        const std::string t  = trim(raw);
        if (t.empty() || t[0] == '!') continue;
        const std::string lo = to_lower(t);

        if (!lo.empty() && (t[0] == '&' || t[0] == '/')) {
            inCell = inAtoms = inSpecies = inKpoints = false; continue;
        }
        if (lo.rfind("cell_parameters", 0) == 0) {
            inCell = true; cellRow = 0; inAtoms = inSpecies = inKpoints = false; continue;
        }
        if (lo.rfind("atomic_positions", 0) == 0) {
            inAtoms = true; inCell = inSpecies = inKpoints = false; continue;
        }
        if (lo.rfind("atomic_species", 0) == 0) {
            inSpecies = true; inCell = inAtoms = inKpoints = false; continue;
        }
        if (lo.rfind("k_points", 0) == 0) {
            inKpoints = true; inCell = inAtoms = inSpecies = false;
            kbuf.str(""); kbuf << t << "\n"; continue;
        }

        if (inCell && cellRow < 3) {
            std::istringstream ss(t); double x,y,z;
            if (ss >> x >> y >> z) qi.cell.row(cellRow++) = Eigen::Vector3d(x,y,z);
            continue;
        }
        if (inAtoms) {
            std::istringstream ss(t); std::string sym; double x,y,z;
            if (ss >> sym >> x >> y >> z)
                qi.atoms.emplace_back(sym, Eigen::Vector3d(x,y,z));
            continue;
        }
        if (inSpecies) { qi.speciesLines.push_back(t); continue; }
        if (inKpoints) { kbuf << t << "\n"; continue; }

        if (lo.find('=') == std::string::npos) continue;
        const auto eq = t.find('=');
        if (lo.find("prefix") != std::string::npos && lo.find("pseudo") == std::string::npos)
            qi.prefix = strip_quotes(stripVal(t, eq));
        else if (lo.find("pseudo_dir") != std::string::npos)
            qi.pseudoDir = strip_quotes(stripVal(t, eq));
        else if (lo.find("outdir") != std::string::npos)
            qi.outdir = strip_quotes(stripVal(t, eq));
        else if (lo.find("ecutwfc") != std::string::npos)
            try { qi.ecutwfc = std::stoi(stripVal(t, eq)); } catch (...) {}
        else if (lo.find("ecutrho") != std::string::npos)
            try { qi.ecutrho = std::stoi(stripVal(t, eq)); } catch (...) {}
        else if (lo.find("occupations") != std::string::npos)
            qi.occupations = strip_quotes(stripVal(t, eq));
        else if (lo.find("smearing") != std::string::npos)
            qi.smearing = strip_quotes(stripVal(t, eq));
        else if (lo.find("degauss") != std::string::npos)
            try { qi.degauss = std::stod(stripVal(t, eq)); } catch (...) {}
        else if (lo.find("ibrav") != std::string::npos && lo.find("celldm") == std::string::npos)
            try { qi.ibrav = std::stoi(stripVal(t, eq)); } catch (...) {}
        else if (lo.find("celldm(") != std::string::npos) {
            const auto lp = lo.find('('), rp = lo.find(')');
            if (lp != std::string::npos && rp != std::string::npos) {
                try {
                    int idx = std::stoi(lo.substr(lp+1, rp-lp-1)) - 1;
                    if (idx >= 0 && idx < 6)
                        qi.celldm[idx] = std::stod(stripVal(t, eq));
                } catch (...) {}
            }
        }
    }
    qi.kpointsBlock = kbuf.str();

    // Build cell from ibrav/celldm when no CELL_PARAMETERS block was found
    if (qi.cell.isZero() && qi.ibrav != 0 && qi.celldm[0] != 0.0) {
        const double a = qi.celldm[0] * bohr2ang;
        switch (qi.ibrav) {
            case 1:  qi.cell << a,0,0,0,a,0,0,0,a; break;  // SC
            case 2:  qi.cell << -a/2,0,a/2, 0,a/2,a/2, -a/2,a/2,0; break;  // FCC
            case 3:  qi.cell << a/2,a/2,a/2, -a/2,a/2,a/2, -a/2,-a/2,a/2; break;  // BCC
            case -3: qi.cell << -a/2,a/2,a/2, a/2,-a/2,a/2, a/2,a/2,-a/2; break;
            case 4: { // Hexagonal
                const double c = qi.celldm[2]*a;
                qi.cell << a,0,0, -a/2,a*std::sqrt(3.0)/2,0, 0,0,c; break;
            }
            case 6: { // Tetragonal P
                const double c = qi.celldm[2]*a;
                qi.cell << a,0,0, 0,a,0, 0,0,c; break;
            }
            case 8: { // Orthorhombic P
                const double b = qi.celldm[1]*a, c = qi.celldm[2]*a;
                qi.cell << a,0,0, 0,b,0, 0,0,c; break;
            }
            default: break;
        }
    }
    const Eigen::Vector3d a1=qi.cell.row(0), a2=qi.cell.row(1), a3=qi.cell.row(2);
    qi.volume = std::abs(a1.dot(a2.cross(a3)));
    return qi;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Crystal family detection: spglib → ibrav heuristic → triclinic
// ─────────────────────────────────────────────────────────────────────────────
static std::string detect_crystal_family(const QEInput& qi, bool verbose = true) {
#ifdef HAVE_SPGLIB
    if (qi.volume > 0.0 && !qi.atoms.empty()) {
        // spglib expects columns as lattice vectors: lat[row][col] = col-th vector's row-th component
        double lat[3][3];
        for (int i=0;i<3;i++) for (int j=0;j<3;j++) lat[i][j] = qi.cell(j,i);
        const int natom = static_cast<int>(qi.atoms.size());
        std::vector<double> pos_flat(natom * 3);
        for (int a=0;a<natom;a++) {
            pos_flat[a*3+0] = qi.atoms[a].second.x();
            pos_flat[a*3+1] = qi.atoms[a].second.y();
            pos_flat[a*3+2] = qi.atoms[a].second.z();
        }
        std::map<std::string,int> symmap; int ncnt = 0;
        std::vector<int> types(natom);
        for (int a=0;a<natom;a++) {
            const auto& sym = qi.atoms[a].first;
            if (!symmap.count(sym)) symmap[sym] = ncnt++;
            types[a] = symmap[sym];
        }
        SpglibDataset* ds = spg_get_dataset(
            lat,
            reinterpret_cast<double(*)[3]>(pos_flat.data()),
            types.data(), natom, 1e-5);
        if (ds) {
            const int sg = ds->spacegroup_number;
            const std::string fam = crystal_family_from_sg(sg);
            if (verbose)
                std::cout << "  spglib: SG=" << sg
                          << " (" << ds->international_symbol << ")"
                          << ", family=" << fam << "\n";
            spg_free_dataset(ds);
            return fam;
        }
    }
#else
    (void)verbose;
#endif
    // Fallback: ibrav heuristic
    const int ib = qi.ibrav;
    if (ib==1||ib==2||ib==3||ib==-3)    return "cubic";
    if (ib==4)                            return "hexagonal";
    if (ib==5||ib==-5)                   return "trigonal";
    if (ib==6||ib==7)                    return "tetragonal";
    if (ib>=8  && ib<=11)                return "orthorhombic";
    if (ib>=12 && ib<=14)                return "monoclinic";
    if (verbose && ib==0)
        std::cerr << "  Warning: ibrav=0 and spglib unavailable — defaulting to triclinic.\n";
    return "triclinic";
}

// ─────────────────────────────────────────────────────────────────────────────
//  Strain patterns
//
//  Convention: Voigt strains (k=0..5) applied with amplitude delta:
//    normal  k=0,1,2: eps_kk = eta[k]*delta
//    shear   k=3,4,5: eps_ij = eta[k]*delta/2  (Voigt gamma_ij = eta[k]*delta)
//
//  Energy relation (exact to 2nd order):
//    A2 = V0/2 * sum_{k,l} C_{kl} * eta_k * eta_l
//
//  => single strain (one non-zero eta_k):   C_{kk} = 2*A2 / (V * eta_k^2)
//  => cross  strain (two non-zero eta_i,j): C_{ij} = (2*A2/V - C_{ii}*eta_i^2 - C_{jj}*eta_j^2)
//                                                     / (2*eta_i*eta_j)
// ─────────────────────────────────────────────────────────────────────────────
static std::vector<StrainPattern> patterns_for_family(const std::string& family) {
    std::vector<StrainPattern> pats;

    auto single = [](int k) -> StrainPattern {
        StrainPattern p; p.name = "e"+std::to_string(k+1);
        p.eta.fill(0.0); p.eta[k] = 1.0; return p;
    };
    auto cross = [](int i, int j) -> StrainPattern {
        StrainPattern p;
        p.name = "e"+std::to_string(i+1)+"e"+std::to_string(j+1);
        p.eta.fill(0.0); p.eta[i] = 1.0; p.eta[j] = 1.0; return p;
    };

    if (family == "cubic") {
        // 3 independent: C11, C12, C44
        pats.push_back(single(0));   // e1  → C11
        pats.push_back(cross(0,1));  // e1e2 → C12 (uses C11=C22)
        pats.push_back(single(3));   // e4  → C44
    }
    else if (family == "hexagonal") {
        // 5 independent: C11=C22, C33, C12, C13=C23, C44=C55; C66=(C11-C12)/2
        pats.push_back(single(0));   // e1  → C11
        pats.push_back(single(2));   // e3  → C33
        pats.push_back(single(3));   // e4  → C44
        pats.push_back(cross(0,1));  // e1e2 → C12
        pats.push_back(cross(0,2));  // e1e3 → C13
    }
    else if (family == "trigonal") {
        // 6 independent: C11=C22, C33, C12, C13=C23, C14=-C24=C56, C44=C55; C66=(C11-C12)/2
        pats.push_back(single(0));
        pats.push_back(single(2));
        pats.push_back(single(3));
        pats.push_back(cross(0,1));
        pats.push_back(cross(0,2));
        pats.push_back(cross(0,3));  // e1e4 → C14
    }
    else if (family == "trigonal_I") {
        // 7 independent: above + C15
        pats.push_back(single(0));
        pats.push_back(single(2));
        pats.push_back(single(3));
        pats.push_back(cross(0,1));
        pats.push_back(cross(0,2));
        pats.push_back(cross(0,3));
        pats.push_back(cross(0,4));  // e1e5 → C15
    }
    else if (family == "tetragonal") {
        // 6 independent: C11=C22, C33, C12, C13=C23, C44=C55, C66
        pats.push_back(single(0));
        pats.push_back(single(2));
        pats.push_back(single(3));
        pats.push_back(single(5));   // e6 → C66
        pats.push_back(cross(0,1));
        pats.push_back(cross(0,2));
    }
    else if (family == "tetragonal_I") {
        // 7 independent: above + C16
        pats.push_back(single(0));
        pats.push_back(single(2));
        pats.push_back(single(3));
        pats.push_back(single(5));
        pats.push_back(cross(0,1));
        pats.push_back(cross(0,2));
        pats.push_back(cross(0,5));  // e1e6 → C16
    }
    else if (family == "orthorhombic") {
        // 9 independent: C11,C22,C33, C12,C13,C23, C44,C55,C66
        for (int k=0;k<6;k++) pats.push_back(single(k));
        pats.push_back(cross(0,1));
        pats.push_back(cross(0,2));
        pats.push_back(cross(1,2));
    }
    else if (family == "monoclinic") {
        // 13 independent (unique axis b): above 9 + C15,C25,C35,C46
        for (int k=0;k<6;k++) pats.push_back(single(k));
        pats.push_back(cross(0,1));
        pats.push_back(cross(0,2));
        pats.push_back(cross(1,2));
        pats.push_back(cross(0,4));  // C15
        pats.push_back(cross(1,4));  // C25
        pats.push_back(cross(2,4));  // C35
        pats.push_back(cross(3,5));  // C46
    }
    else {
        // triclinic: 21 independent — 6 single + 15 cross
        for (int k=0;k<6;k++) pats.push_back(single(k));
        for (int i=0;i<6;i++)
            for (int j=i+1;j<6;j++)
                pats.push_back(cross(i,j));
    }
    return pats;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Apply Voigt strain to cell (row vectors, Å)
// ─────────────────────────────────────────────────────────────────────────────
static Eigen::Matrix3d apply_strain(const Eigen::Matrix3d& cell,
                                    const std::array<double,6>& eta,
                                    double delta) {
    Eigen::Matrix3d eps = Eigen::Matrix3d::Zero();
    eps(0,0)=eta[0]*delta; eps(1,1)=eta[1]*delta; eps(2,2)=eta[2]*delta;
    eps(1,2)=eps(2,1)=0.5*eta[3]*delta;
    eps(0,2)=eps(2,0)=0.5*eta[4]*delta;
    eps(0,1)=eps(1,0)=0.5*eta[5]*delta;
    return ((Eigen::Matrix3d::Identity()+eps)*cell.transpose()).transpose();
}

// ─────────────────────────────────────────────────────────────────────────────
//  Write deformed QE SCF input (always ibrav=0 + CELL_PARAMETERS angstrom)
// ─────────────────────────────────────────────────────────────────────────────
static void write_deformed_input(const std::string& outPath,
                                  const QEInput& base,
                                  const Eigen::Matrix3d& newCell,
                                  const std::string& outdir,
                                  const std::string& overridePseudo = "") {
    const std::string& pd = overridePseudo.empty() ? base.pseudoDir : overridePseudo;
    std::ofstream out(outPath);
    if (!out.is_open()) throw std::runtime_error("Cannot open: " + outPath);
    out << "&CONTROL\n"
        << "  calculation = 'scf',\n"
        << "  prefix = '" << base.prefix << "',\n"
        << "  outdir = '" << outdir << "',\n"
        << "  pseudo_dir = '" << pd << "',\n"
        << "  verbosity = 'low'\n"
        << "/\n\n"
        << "&SYSTEM\n"
        << "  ibrav = 0,\n"
        << "  nat = " << base.atoms.size() << ",\n"
        << "  ntyp = " << base.speciesLines.size() << ",\n"
        << "  ecutwfc = " << base.ecutwfc << ",\n"
        << "  ecutrho = " << base.ecutrho << "\n";
    if (!base.occupations.empty())
        out << "  occupations = '" << base.occupations << "',\n";
    if (!base.smearing.empty())
        out << "  smearing = '" << base.smearing << "',\n";
    if (base.degauss > 0.0)
        out << "  degauss = " << base.degauss << ",\n";
    out << "/\n\n"
        << "&ELECTRONS\n"
        << "  conv_thr = 1.0d-10\n"
        << "/\n\n"
        << "ATOMIC_SPECIES\n";
    for (const auto& s : base.speciesLines) out << "  " << s << "\n";
    out << "\nCELL_PARAMETERS angstrom\n" << std::fixed << std::setprecision(10);
    for (int r=0;r<3;r++)
        out << "  " << std::setw(16) << newCell(r,0)
            << "  " << std::setw(16) << newCell(r,1)
            << "  " << std::setw(16) << newCell(r,2) << "\n";
    out << "\nATOMIC_POSITIONS crystal\n";
    for (const auto& [sym, frac] : base.atoms)
        out << "  " << sym
            << "  " << std::setw(14) << frac.x()
            << "  " << std::setw(14) << frac.y()
            << "  " << std::setw(14) << frac.z() << "\n";
    out << "\n" << base.kpointsBlock;
}

// ─────────────────────────────────────────────────────────────────────────────
//  PRE-PROCESSING: generate_elastic_inputs
// ─────────────────────────────────────────────────────────────────────────────
void generate_elastic_inputs(const std::string& scfTemplatePath,
                              const std::string& outDir,
                              int nDeltas, double maxDelta) {
    if (nDeltas < 3 || nDeltas % 2 == 0)
        throw std::runtime_error("nDeltas must be odd and >= 3.");
    if (maxDelta <= 0.0 || maxDelta >= 0.2)
        throw std::runtime_error("maxDelta must be in (0, 0.2).");

    const QEInput base = parse_qe_input(scfTemplatePath);
    if (base.volume <= 0.0)
        throw std::runtime_error(
            "Could not determine unit cell volume from: " + scfTemplatePath +
            "\nEnsure CELL_PARAMETERS angstrom or ibrav + celldm(1) are present.");

    std::cout << "Detecting crystal family...\n";
    const std::string family  = detect_crystal_family(base, /*verbose=*/true);
    const auto        patterns = patterns_for_family(family);

    std::vector<double> deltas;
    for (int i=0;i<nDeltas;i++) {
        double d = -maxDelta + i*(2.0*maxDelta/(nDeltas-1));
        if (std::abs(d) < 1e-12) d = 0.0;
        deltas.push_back(d);
    }

    // Resolve pseudo_dir to absolute path
    std::string resolvedPseudo = base.pseudoDir;
    if (!resolvedPseudo.empty() && resolvedPseudo[0] != '/') {
        const fs::path tmplDir = fs::path(scfTemplatePath).parent_path();
        resolvedPseudo = fs::weakly_canonical(tmplDir / resolvedPseudo).string();
    }

    fs::create_directories(outDir);
    {
        std::ofstream meta(outDir + "/elastic_setup.dat");
        meta << "ndeltas "  << nDeltas << "\n"
             << "max_delta " << std::fixed << std::setprecision(6) << maxDelta << "\n";
    }

    int total = 0;
    for (const auto& pat : patterns) {
        fs::create_directories(outDir + "/" + pat.name);
        for (double d : deltas) {
            std::ostringstream dn;
            dn << std::fixed << std::setprecision(4) << (d>=0?"p":"m") << std::abs(d);
            const std::string dDir = outDir + "/" + pat.name + "/" + dn.str();
            fs::create_directories(dDir);
            write_deformed_input(dDir + "/" + base.prefix + ".in", base,
                                 apply_strain(base.cell, pat.eta, d),
                                 "./tmp", resolvedPseudo);
            ++total;
        }
    }

    std::cout << "Crystal family   : " << family << "\n"
              << "Volume (Ang^3)   : " << std::fixed << std::setprecision(4) << base.volume << "\n"
              << "Patterns         : " << patterns.size() << "\n"
              << "Strains/pattern  : " << nDeltas << " (" << -maxDelta << " to +" << maxDelta << ")\n"
              << "Total inputs     : " << total << "\n"
              << "Output directory : " << outDir << "\n\n"
              << "Run all calculations with:\n"
              << "  for d in " << outDir << "/*/*/; do\n"
              << "    (cd \"$d\" && pw.x < " << base.prefix << ".in > "
                                             << base.prefix << ".out)\n"
              << "  done\n";
}

// ─────────────────────────────────────────────────────────────────────────────
//  POST-PROCESSING helpers
// ─────────────────────────────────────────────────────────────────────────────

static double fit_quadratic_a2(const std::vector<double>& x,
                                const std::vector<double>& y) {
    const int n = static_cast<int>(x.size());
    if (n < 3) throw std::runtime_error("Need >= 3 points for quadratic fit.");
    Eigen::MatrixXd A(n,3); Eigen::VectorXd b(n);
    for (int i=0;i<n;i++) { A(i,0)=1.0; A(i,1)=x[i]; A(i,2)=x[i]*x[i]; b[i]=y[i]; }
    return A.colPivHouseholderQr().solve(b)[2];
}

// ─────────────────────────────────────────────────────────────────────────────
//  Assemble 6x6 C (in Pa) from A2 map using the universal two-pass algorithm:
//    Pass 1  — single strains  → diagonal C_kk
//    Pass 1b — propagate symmetry-equivalent diagonal entries so that cross-
//              strain formulas can use the correct C_jj (not zero)
//    Pass 2  — cross  strains  → off-diagonal C_ij
// ─────────────────────────────────────────────────────────────────────────────
static Eigen::Matrix<double,6,6>
assemble_C(const std::vector<StrainPattern>& patterns,
           const std::map<std::string,double>& a2map,
           double V_m3,
           const std::string& family) {
    Eigen::Matrix<double,6,6> C = Eigen::Matrix<double,6,6>::Zero();

    // Pass 1: single-strain patterns → diagonal
    for (const auto& pat : patterns) {
        if (!a2map.count(pat.name)) continue;
        std::vector<int> nz;
        for (int k=0;k<6;k++) if (std::abs(pat.eta[k])>1e-10) nz.push_back(k);
        if (nz.size() == 1) {
            const int k = nz[0];
            C(k,k) = 2.0 * a2map.at(pat.name)*eV2J / (V_m3 * pat.eta[k]*pat.eta[k]);
        }
    }

    // Pass 1b: propagate symmetry-equivalent diagonal entries.
    // When only one representative pattern is computed per symmetry class,
    // the other equivalent entries stay zero without this step, causing
    // incorrect cross-term computation in pass 2.
    auto prop2 = [&](int a, int b) {
        if (C(a,a) != 0.0 && C(b,b) == 0.0) C(b,b) = C(a,a);
        else if (C(b,b) != 0.0 && C(a,a) == 0.0) C(a,a) = C(b,b);
    };
    auto prop3 = [&](int a, int b, int c) {
        double v = (C(a,a) != 0.0) ? C(a,a) : ((C(b,b) != 0.0) ? C(b,b) : C(c,c));
        if (v != 0.0) { C(a,a) = C(b,b) = C(c,c) = v; }
    };
    if (family == "cubic") {
        prop3(0, 1, 2);  // C11=C22=C33
        prop3(3, 4, 5);  // C44=C55=C66
    } else if (family == "hexagonal" || family == "trigonal" ||
               family == "trigonal_I" || family == "tetragonal" ||
               family == "tetragonal_I") {
        prop2(0, 1);  // C11=C22
        prop2(3, 4);  // C44=C55
    }

    // Pass 2: cross-strain patterns → off-diagonal
    for (const auto& pat : patterns) {
        if (!a2map.count(pat.name)) continue;
        std::vector<int> nz;
        for (int k=0;k<6;k++) if (std::abs(pat.eta[k])>1e-10) nz.push_back(k);
        if (nz.size() == 2) {
            const int i=nz[0], j=nz[1];
            const double a2J = a2map.at(pat.name)*eV2J;
            const double rhs = 2.0*a2J/V_m3
                               - C(i,i)*pat.eta[i]*pat.eta[i]
                               - C(j,j)*pat.eta[j]*pat.eta[j];
            C(i,j) = C(j,i) = rhs / (2.0*pat.eta[i]*pat.eta[j]);
        }
    }
    return C;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Apply crystal-symmetry constraints (averaging symmetry-equivalent entries)
// ─────────────────────────────────────────────────────────────────────────────
static void apply_symmetry(Eigen::Matrix<double,6,6>& C, const std::string& family) {
    auto avg2=[&](int i,int j,int k,int l){
        double v=0.5*(C(i,j)+C(k,l)); C(i,j)=C(j,i)=C(k,l)=C(l,k)=v; };
    auto avg3=[&](int i,int j,int k,int l,int m,int n){
        double v=(C(i,j)+C(k,l)+C(m,n))/3.0;
        C(i,j)=C(j,i)=C(k,l)=C(l,k)=C(m,n)=C(n,m)=v; };

    if (family == "cubic") {
        // Propagate off-diagonal: only C(0,1) was computed from e1e2;
        // set C13=C23=C12 before averaging.
        {
            double v = (C(0,1) != 0.0) ? C(0,1) : ((C(0,2) != 0.0) ? C(0,2) : C(1,2));
            if (v != 0.0) { C(0,1)=C(1,0)=C(0,2)=C(2,0)=C(1,2)=C(2,1)=v; }
        }
        avg3(0,0,1,1,2,2); avg3(0,1,0,2,1,2); avg3(3,3,4,4,5,5);
        for (int i=0;i<6;i++) for (int j=0;j<6;j++) {
            bool keep=(i==j&&i<=2)||(i!=j&&i<=2&&j<=2)||(i==j&&i>=3);
            if (!keep) C(i,j)=0.0;
        }
    }
    else if (family == "hexagonal") {
        // Propagate C23=C13 before averaging
        if (C(0,2) != 0.0 && C(1,2) == 0.0) C(1,2) = C(2,1) = C(0,2);
        else if (C(1,2) != 0.0 && C(0,2) == 0.0) C(0,2) = C(2,0) = C(1,2);
        avg2(0,0,1,1); avg2(0,2,1,2); avg2(3,3,4,4);
        C(5,5) = 0.5*(C(0,0)-C(0,1));
        // Zero off-axis coupling
        for (int i=0;i<3;i++) for (int j=3;j<6;j++) C(i,j)=C(j,i)=0.0;
        for (int i=3;i<6;i++) for (int j=3;j<6;j++) if (i!=j) C(i,j)=0.0;
    }
    else if (family == "tetragonal" || family == "tetragonal_I") {
        // Propagate C23=C13 before averaging
        if (C(0,2) != 0.0 && C(1,2) == 0.0) C(1,2) = C(2,1) = C(0,2);
        else if (C(1,2) != 0.0 && C(0,2) == 0.0) C(0,2) = C(2,0) = C(1,2);
        avg2(0,0,1,1); avg2(0,2,1,2); avg2(3,3,4,4);
        if (family == "tetragonal_I") {
            double v=0.5*(std::abs(C(0,5))+std::abs(C(1,5)));
            C(0,5)=C(5,0)=v; C(1,5)=C(5,1)=-v;
        } else {
            C(0,5)=C(5,0)=C(1,5)=C(5,1)=0.0;
        }
    }
    else if (family == "trigonal" || family == "trigonal_I") {
        // Propagate C23=C13 before averaging
        if (C(0,2) != 0.0 && C(1,2) == 0.0) C(1,2) = C(2,1) = C(0,2);
        else if (C(1,2) != 0.0 && C(0,2) == 0.0) C(0,2) = C(2,0) = C(1,2);
        avg2(0,0,1,1); avg2(0,2,1,2); avg2(3,3,4,4);
        C(5,5) = 0.5*(C(0,0)-C(0,1));
        // C14=-C24=C56
        double c14=(C(0,3)-C(1,3)+C(4,5))/3.0;
        C(0,3)=C(3,0)=c14; C(1,3)=C(3,1)=-c14; C(4,5)=C(5,4)=c14;
    }
    // Enforce symmetry numerically
    for (int i=0;i<6;i++) for (int j=0;j<6;j++) C(i,j)=C(j,i)=0.5*(C(i,j)+C(j,i));
}

// ─────────────────────────────────────────────────────────────────────────────
//  Mechanical stability — Born criteria per crystal family + eigenvalue check
// ─────────────────────────────────────────────────────────────────────────────
static void check_stability(const Eigen::Matrix<double,6,6>& C,
                             const std::string& family,
                             bool& stable,
                             std::vector<std::string>& notes) {
    stable = true; notes.clear();

    // General: all eigenvalues of C must be positive
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,6,6>> es(C);
    const double minEig = es.eigenvalues().minCoeff();
    if (minEig < 0.0) {
        stable = false;
        notes.push_back("Eigenvalue check FAILED: min eigenvalue = "
                        + std::to_string(minEig) + " GPa");
    }

    const double C11=C(0,0),C22=C(1,1),C33=C(2,2);
    const double C12=C(0,1),C13=C(0,2),C23=C(1,2);
    const double C44=C(3,3),C55=C(4,4),C66=C(5,5);

    auto flag=[&](bool cond,const std::string& msg){
        if(!cond){stable=false;notes.push_back("FAILED: "+msg);}
    };

    if (family == "cubic") {
        flag(C11-std::abs(C12)>0,  "Cubic: C11 - |C12| > 0");
        flag(C11+2*C12>0,           "Cubic: C11 + 2*C12 > 0");
        flag(C44>0,                  "Cubic: C44 > 0");
    }
    else if (family == "hexagonal") {
        flag(C11>std::abs(C12),          "Hexagonal: C11 > |C12|");
        flag(C33*(C11+C12)>2*C13*C13,    "Hexagonal: C33*(C11+C12) > 2*C13^2");
        flag(C44>0,                       "Hexagonal: C44 > 0");
    }
    else if (family == "trigonal" || family == "trigonal_I") {
        flag(C11>std::abs(C12),          "Trigonal: C11 > |C12|");
        flag(C33*(C11+C12)>2*C13*C13,    "Trigonal: C33*(C11+C12) > 2*C13^2");
        flag(C44>0,                       "Trigonal: C44 > 0");
        if (C44>1e-10)
            flag((C11-C12)*C44>2*C(0,3)*C(0,3),
                 "Trigonal: (C11-C12)*C44 > 2*C14^2");
    }
    else if (family == "tetragonal" || family == "tetragonal_I") {
        flag(C11>std::abs(C12),          "Tetragonal: C11 > |C12|");
        flag(2*C13*C13<C33*(C11+C12),    "Tetragonal: 2*C13^2 < C33*(C11+C12)");
        flag(C44>0,                       "Tetragonal: C44 > 0");
        flag(C66>0,                       "Tetragonal: C66 > 0");
    }
    else if (family == "orthorhombic") {
        flag(C11>0,"Orthorhombic: C11>0"); flag(C22>0,"Orthorhombic: C22>0");
        flag(C33>0,"Orthorhombic: C33>0"); flag(C44>0,"Orthorhombic: C44>0");
        flag(C55>0,"Orthorhombic: C55>0"); flag(C66>0,"Orthorhombic: C66>0");
        flag(C11*C22>C12*C12, "Orthorhombic: C11*C22 > C12^2");
        const double det3=C11*(C22*C33-C23*C23)-C12*(C12*C33-C23*C13)+C13*(C12*C23-C22*C13);
        flag(det3>0, "Orthorhombic: det(C_3x3) > 0");
    }
    else if (family == "monoclinic") {
        flag(C11>0,"Mono: C11>0"); flag(C22>0,"Mono: C22>0"); flag(C33>0,"Mono: C33>0");
        flag(C44>0,"Mono: C44>0"); flag(C55>0,"Mono: C55>0"); flag(C66>0,"Mono: C66>0");
        flag(C11*C22>C12*C12, "Mono: C11*C22 > C12^2");
        const double det3=C11*(C22*C33-C23*C23)-C12*(C12*C33-C23*C13)+C13*(C12*C23-C22*C13);
        flag(det3>0, "Mono: det(C_3x3) > 0");
    }
    // triclinic: eigenvalue check is the sole criterion (already done above)
}

// ─────────────────────────────────────────────────────────────────────────────
//  POST-PROCESSING: compute_elastic_properties
// ─────────────────────────────────────────────────────────────────────────────
ElasticResults compute_elastic_properties(const std::string& scfTemplatePath,
                                          const std::string& outDir) {
    // --- Read setup file (only the user-chosen run parameters) ---
    int ndeltas=0; double maxDelta=0.0;
    {
        std::ifstream meta(outDir+"/elastic_setup.dat");
        if (!meta.is_open())
            throw std::runtime_error(
                "Cannot open "+outDir+"/elastic_setup.dat\n"
                "Run 'qepp elastic ...' first.");
        std::string line;
        while (std::getline(meta,line)) {
            std::istringstream ss(line); std::string key; ss>>key;
            if (key=="ndeltas")   ss>>ndeltas;
            if (key=="max_delta") ss>>maxDelta;
        }
    }

    // --- Derive structure-dependent quantities from the template ---
    const QEInput base = parse_qe_input(scfTemplatePath);
    const std::string family = detect_crystal_family(base, /*verbose=*/false);
    const std::vector<StrainPattern> patterns = patterns_for_family(family);
    const double volume_ang3 = (base.volume > 0.0) ? base.volume
        : throw std::runtime_error("Volume is zero — check cell parameters in template.");
    const double V_m3 = volume_ang3 * Ang2m*Ang2m*Ang2m;

    std::vector<double> deltas;
    for (int i=0;i<ndeltas;i++) {
        double d = -maxDelta+i*(2.0*maxDelta/(ndeltas-1));
        if (std::abs(d)<1e-12) d=0.0;
        deltas.push_back(d);
    }

    // --- Collect energies and fit A2 per pattern ---
    std::map<std::string,double> a2map;
    for (const auto& pat : patterns) {
        std::vector<double> dx, dy;
        for (double d : deltas) {
            std::ostringstream dn;
            dn << std::fixed << std::setprecision(4) << (d>=0?"p":"m") << std::abs(d);
            const std::string outFile = outDir+"/"+pat.name+"/"+dn.str()+"/"+base.prefix+".out";
            try { dx.push_back(d); dy.push_back(extract_total_energy(outFile)); }
            catch (const std::exception& ex) {
                std::cerr << "  Warning: skipping " << outFile << " — " << ex.what() << "\n";
            }
        }
        if (static_cast<int>(dx.size()) < 3)
            throw std::runtime_error(
                "Not enough converged data for pattern '"+pat.name+"' (found "
                +std::to_string(dx.size())+", need >=3).");
        a2map[pat.name] = fit_quadratic_a2(dx, dy);
    }

    // --- Assemble and symmetrize C (GPa) ---
    Eigen::Matrix<double,6,6> C_Pa = assemble_C(patterns, a2map, V_m3, family);
    apply_symmetry(C_Pa, family);

    ElasticResults res;
    res.crystalFamily = family;
    res.C = C_Pa / 1.0e9;
    res.S = res.C.inverse();

    // --- VRH bulk and shear moduli ---
    const auto& C = res.C; const auto& S = res.S;
    res.KV = (C(0,0)+C(1,1)+C(2,2)+2*(C(0,1)+C(0,2)+C(1,2)))/9.0;
    res.GV = (C(0,0)+C(1,1)+C(2,2)-C(0,1)-C(0,2)-C(1,2)+3*(C(3,3)+C(4,4)+C(5,5)))/15.0;
    const double KRi = S(0,0)+S(1,1)+S(2,2)+2*(S(0,1)+S(0,2)+S(1,2));
    res.KR = (std::abs(KRi)>1e-30)?1.0/KRi:0.0;
    const double GRi = 4*(S(0,0)+S(1,1)+S(2,2))-4*(S(0,1)+S(0,2)+S(1,2))+3*(S(3,3)+S(4,4)+S(5,5));
    res.GR = (std::abs(GRi)>1e-30)?15.0/GRi:0.0;
    res.KH=0.5*(res.KV+res.KR); res.GH=0.5*(res.GV+res.GR);
    res.EH=9*res.KH*res.GH/(3*res.KH+res.GH);
    res.nuH=(3*res.KH-2*res.GH)/(2*(3*res.KH+res.GH));
    const double d=C(0,0)-C(0,1);
    res.AZ=(std::abs(d)>1e-10)?2*C(3,3)/d:1.0;

    // --- Sound velocities ---
    double massAmu=0.0;
    for (const auto& [sym,frac]:base.atoms) massAmu+=atomic_mass(sym);
    const double rho=massAmu*amuToKg/V_m3;
    const double KH_Pa=res.KH*GPa2Pa, GH_Pa=res.GH*GPa2Pa;
    res.vL=std::sqrt((KH_Pa+4.0/3.0*GH_Pa)/rho);
    res.vT=std::sqrt(GH_Pa/rho);
    res.vM=std::pow((1.0/3.0)*(2.0/std::pow(res.vT,3)+1.0/std::pow(res.vL,3)),-1.0/3.0);
    res.vDebye=res.vM;

    // --- Debye temperature ---
    const int nA=static_cast<int>(base.atoms.size());
    res.thetaDebye=(hbar/kB)*std::pow(6*M_PI*M_PI*nA/V_m3,1.0/3.0)*res.vM;

    // --- Lame constants (Hill) ---
    res.lambdaH = res.KH - 2.0/3.0*res.GH;
    res.muH     = res.GH;

    // --- Mechanical character ---
    res.pughRatio      = (res.GH > 1e-10) ? res.KH/res.GH : 0.0;
    res.cauchyPressure = res.C(0,1) - res.C(3,3);

    // --- Anisotropy indices ---
    res.AU = (res.GR>1e-10 && res.KR>1e-10) ? 5.0*res.GV/res.GR + res.KV/res.KR - 6.0 : 0.0;
    res.AB = (res.KV+res.KR>1e-10) ? 100.0*(res.KV-res.KR)/(res.KV+res.KR) : 0.0;
    res.AG = (res.GV+res.GR>1e-10) ? 100.0*(res.GV-res.GR)/(res.GV+res.GR) : 0.0;

    // --- Hardness estimates ---
    {
        const double k = (res.KH > 1e-10) ? res.GH/res.KH : 0.0;
        const double k2G = k*k*res.GH;
        res.hardnessChen = (k2G > 0.0) ? 2.0*std::pow(k2G, 0.585) - 3.0 : 0.0;
        res.hardnessTian = (k > 0.0 && res.GH > 0.0)
                         ? 0.92*std::pow(k, 1.137)*std::pow(res.GH, 0.708) : 0.0;
    }

    // --- Grüneisen parameter (acoustic, Poisson's ratio route) ---
    {
        const double denom = 2.0 - 3.0*res.nuH;
        res.gruneisen = (std::abs(denom) > 1e-10)
                      ? 3.0*(1.0 + res.nuH)/(2.0*denom) : 0.0;
    }

    // --- Minimum thermal conductivity (Cahill-Watson-Pohl model) ---
    {
        const double n_at = nA / V_m3;   // atoms/m³
        res.kMinClarke = std::pow(M_PI/6.0, 1.0/3.0)
                       * kB * std::pow(n_at, 2.0/3.0)
                       * (2.0*res.vT + res.vL);
    }

    // --- Empirical melting temperature (Fine, Brown & Marcus, Scr. Metall. 1984) ---
    res.meltingTemp = 553.0 + 5.91*res.C(0,0);  // C11 in GPa → Tm in K

    // --- Additional hardness models ---
    res.hardnessTeter = 0.151 * res.GV;  // MRS Bull. 23 (1998) 22
    res.hardnessNiu   = (res.nuH < 0.5 && res.nuH > -1.0)
                      ? res.GH * (1.0 - 2.0*res.nuH) / 3.0 : 0.0;

    // --- Directional Young's moduli from compliance (Nye 1957) ---
    {
        res.E_x = (std::abs(res.S(0,0)) > 1e-15) ? 1.0/res.S(0,0) : 0.0;
        res.E_y = (std::abs(res.S(1,1)) > 1e-15) ? 1.0/res.S(1,1) : 0.0;
        res.E_z = (std::abs(res.S(2,2)) > 1e-15) ? 1.0/res.S(2,2) : 0.0;
        if (family == "cubic") {
            const double s11 = res.S(0,0), s12 = res.S(0,1), s44 = res.S(3,3);
            const double inv110 = 0.5*s11 + 0.5*s12 + 0.25*s44;
            const double inv111 = (s11 + 2.0*s12 + s44) / 3.0;
            res.E110 = (std::abs(inv110) > 1e-15) ? 1.0/inv110 : 0.0;
            res.E111 = (std::abs(inv111) > 1e-15) ? 1.0/inv111 : 0.0;
        }
    }

    // --- Linear compressibility (1/GPa) under hydrostatic pressure ---
    res.beta_x = res.S(0,0) + res.S(0,1) + res.S(0,2);
    res.beta_y = res.S(1,0) + res.S(1,1) + res.S(1,2);
    res.beta_z = res.S(2,0) + res.S(2,1) + res.S(2,2);

    // --- Per-plane shear anisotropy (Chung & Buessem, J. Appl. Phys. 1967) ---
    {
        auto sdiv = [](double a, double b) -> double {
            return (std::abs(b) > 1e-10) ? a/b : 0.0;
        };
        res.A1_plane = sdiv(4.0*res.C(3,3), res.C(0,0)+res.C(2,2)-2.0*res.C(0,2));
        res.A2_plane = sdiv(4.0*res.C(4,4), res.C(1,1)+res.C(2,2)-2.0*res.C(1,2));
        res.A3_plane = sdiv(4.0*res.C(5,5), res.C(0,0)+res.C(1,1)-2.0*res.C(0,1));
    }

    // --- Slack lattice thermal conductivity at 300 K (Slack, JPCS 34 (1973) 321) ---
    {
        const double T_slack   = 300.0;
        const double M_avg_amu = massAmu / nA;      // average mass per atom (amu)
        const double Va_A3     = volume_ang3 / nA;  // volume per atom (Å³)
        const double gamma2    = res.gruneisen * res.gruneisen;
        if (gamma2 > 1e-10 && res.thetaDebye > 1.0) {
            res.kSlack = 3.1e-6 * M_avg_amu
                       * std::pow(res.thetaDebye, 3.0)
                       * std::cbrt(Va_A3)
                       / (gamma2 * std::pow(static_cast<double>(nA), 2.0/3.0) * T_slack);
        }
    }

    // --- Stability ---
    check_stability(res.C, family, res.mechanicallyStable, res.stabilityNotes);
    return res;
}

}  // namespace qe
