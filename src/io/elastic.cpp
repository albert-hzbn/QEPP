#include "qe/elastic.hpp"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

#include "qe/utils.hpp"

namespace qe {

namespace {

std::string shell_quote(const std::string& text) {
   std::string out = "'";
   for (char c : text) {
      if (c == '\'') out += "'\\''";
      else out.push_back(c);
   }
   out += "'";
   return out;
}

bool contains_job_done(const std::filesystem::path& path) {
   std::ifstream in(path);
   if (!in.is_open()) return false;
   std::string line;
   while (std::getline(in, line)) {
      if (line.find("JOB DONE") != std::string::npos) return true;
   }
   return false;
}

}  // namespace

static constexpr double kRyToEv = 13.605693122;

double extract_total_energy(const std::string& path) {
   std::ifstream in(path);
   if (!in.is_open()) throw std::runtime_error("Cannot open QE output: " + path);
   double energy = std::numeric_limits<double>::quiet_NaN();
   std::string line;
   while (std::getline(in, line)) {
      const std::string lo = to_lower(trim(line));
      if (lo.find("!    total energy") != std::string::npos ||
         lo.find("!  total energy")  != std::string::npos) {
         const auto eq = line.find('=');
         if (eq != std::string::npos) {
            std::istringstream ss(line.substr(eq + 1));
            double val;
            std::string unit;
            if (ss >> val >> unit)
               energy = (to_lower(unit) == "ry") ? val * kRyToEv : val;
         }
      }
   }
   if (std::isnan(energy))
      throw std::runtime_error("Total energy not found in: " + path);
   return energy;
}

   void run_elastic_dataset(const std::string& scfTemplatePath,
                      const std::string& outDir,
                      const QeParallelOptions& parallel) {
      namespace fs = std::filesystem;

      const fs::path tmpl(scfTemplatePath);
      const fs::path root(outDir);
      if (!fs::exists(tmpl))
         throw std::runtime_error("SCF template not found: " + scfTemplatePath);
      if (!fs::is_directory(root))
         throw std::runtime_error("Elastic directory not found: " + outDir);

      const std::string launch = qe_parallel_args(parallel);
      const fs::path eqOut = tmpl.parent_path() / "qe.out";
      if (!contains_job_done(eqOut)) {
         std::cout << "Equilibrium SCF: running...\n";
         const std::string cmd = "cd " + shell_quote(tmpl.parent_path().string()) +
            " && " + launch + " pw.x -input " + shell_quote(tmpl.filename().string()) +
            " > qe.out 2>&1";
         const int rc = std::system(cmd.c_str());
         (void)rc;
         if (!contains_job_done(eqOut))
            throw std::runtime_error("Equilibrium SCF failed: " + eqOut.string());
         std::cout << "Equilibrium SCF: done\n";
      } else {
         std::cout << "Equilibrium SCF: done\n";
      }

      for (const auto& patt : fs::directory_iterator(root)) {
         if (!patt.is_directory()) continue;
         for (const auto& strain : fs::directory_iterator(patt.path())) {
            if (!strain.is_directory()) continue;

            fs::path infile;
            for (const auto& file : fs::directory_iterator(strain.path())) {
               if (file.is_regular_file() && file.path().extension() == ".in") {
                  infile = file.path();
                  break;
               }
            }
            if (infile.empty()) continue;

            const fs::path outfile = strain.path() / (infile.stem().string() + ".out");
            if (contains_job_done(outfile)) continue;

            std::cout << "Running " << patt.path().filename().string()
                    << "/" << strain.path().filename().string() << "...\n";
            const std::string cmd = "cd " + shell_quote(strain.path().string()) +
               " && " + launch + " pw.x < " + shell_quote(infile.filename().string()) +
               " > " + shell_quote(outfile.filename().string()) + " 2>&1";
            const int rc = std::system(cmd.c_str());
            (void)rc;
            if (!contains_job_done(outfile))
               throw std::runtime_error("Elastic strain SCF failed: " + outfile.string());
         }
      }
   }

void print_elastic_report(const ElasticResults& res,
                          const std::string& crystalFamily,
                          const std::string& outFile) {
    std::ofstream fout;
    if (!outFile.empty()) {
        fout.open(outFile);
        if (fout.is_open())
            std::cout << "Saving report to: " << outFile << "\n";
        else
            std::cerr << "Warning: cannot open '" << outFile << "'\n";
    }

    auto W  = [&](const std::string& s){ std::cout<<s; if(fout.is_open()) fout<<s; };
    auto sep= [&](int n=72){ W(std::string(n,'-')+"\n"); };
    std::ostringstream buf;

    // Header
    sep();
    W("  ELASTIC CONSTANTS REPORT\n");
    W("  Crystal family : " + crystalFamily + "\n");
    sep();

    // Stiffness matrix
    W("\nStiffness matrix C_ij (GPa):\n");
    buf.str(""); buf<<std::fixed<<std::setprecision(3);
    for (int i=0;i<6;i++) {
        buf<<"  ";
        for (int j=0;j<6;j++) buf<<std::setw(10)<<res.C(i,j);
        buf<<"\n";
    }
    W(buf.str());

    // Compliance matrix
    W("\nCompliance matrix S_ij (1/GPa):\n");
    buf.str(""); buf<<std::setprecision(6);
    for (int i=0;i<6;i++) {
        buf<<"  ";
        for (int j=0;j<6;j++) buf<<std::setw(12)<<res.S(i,j);
        buf<<"\n";
    }
    W(buf.str());

    // Elastic moduli (VRH)
    sep();
    W("  Elastic moduli - Voigt / Reuss / Hill average (GPa)\n");
    W("  Ref: Voigt, Lehrbuch der Kristallphysik (1928);\n"
      "       Reuss, ZAMM 9 (1929) 49;\n"
      "       Hill, Proc. Phys. Soc. A 65 (1952) 349\n");
    sep();
    buf.str(""); buf<<std::fixed<<std::setprecision(3);
    buf<<"  Bulk modulus    K_V="<<std::setw(8)<<res.KV
       <<"  K_R="<<std::setw(8)<<res.KR
       <<"  K_H="<<std::setw(8)<<res.KH<<"\n"
       <<"  Shear modulus   G_V="<<std::setw(8)<<res.GV
       <<"  G_R="<<std::setw(8)<<res.GR
       <<"  G_H="<<std::setw(8)<<res.GH<<"\n"
       <<"  Young's modulus E_H="<<std::setw(8)<<res.EH<<" GPa\n"
       <<"  Poisson's ratio v_H="<<std::setw(8)<<res.nuH<<"\n"
       <<"  Zener anisotropy AZ="<<std::setw(8)<<res.AZ;
    if (std::abs(res.AZ-1.0)<0.05) buf<<"  (near-isotropic)";
    buf<<"  [Zener, Elasticity & Anelasticity of Metals (1948)]\n";
    W(buf.str());

    // Lame constants
    sep();
    W("  Lame constants (Hill, GPa)\n");
    W("  Ref: Landau & Lifshitz, Theory of Elasticity, 3rd ed. (1986)\n");
    sep();
    buf.str(""); buf<<std::fixed<<std::setprecision(3);
    buf<<"  Lambda (1st Lame, lambda = K - 2G/3) = "<<std::setw(10)<<res.lambdaH<<" GPa\n"
       <<"  Mu    (shear,    mu = G)            = "<<std::setw(10)<<res.muH    <<" GPa\n";
    W(buf.str());

    // Mechanical character & anisotropy
    sep();
    W("  Mechanical character & anisotropy indices\n");
    sep();
    buf.str(""); buf<<std::fixed<<std::setprecision(3);
    buf<<"  Pugh's ratio K/G = "<<std::setw(8)<<res.pughRatio
       <<(res.pughRatio > 1.75 ? "  (ductile)" : "  (brittle)")
       <<"  [Pugh, Phil. Mag. 45 (1954) 823]\n";
    buf<<"  Cauchy pressure C12-C44 = "<<std::setw(8)<<res.cauchyPressure<<" GPa"
       <<(res.cauchyPressure > 0 ? "  (metallic/ionic)" : "  (covalent)")
       <<"  [Pettifor, Mater. Sci. Tech. 8 (1992) 345]\n";
    buf<<std::setprecision(4)
       <<"  Universal anisotropy A^U = "<<std::setw(9)<<res.AU
       <<"  (0 = isotropic)  [Ranganathan & Ostoja-Starzewski, PRL 101 (2008) 055504]\n";
    buf<<std::setprecision(2)
       <<"  Percent anisotropy A_B = "<<std::setw(8)<<res.AB<<" %\n"
       <<"  Percent anisotropy A_G = "<<std::setw(8)<<res.AG<<" %\n";
    buf<<"  Per-plane shear anisotropy (A=1 isotropic):"
       <<"  [Chung & Buessem, J. Appl. Phys. 38 (1967) 2535]\n";
    buf<<std::setprecision(4)
       <<"    A1 {100} = "<<std::setw(8)<<res.A1_plane
       <<"    A2 {010} = "<<std::setw(8)<<res.A2_plane
       <<"    A3 {001} = "<<std::setw(8)<<res.A3_plane<<"\n";
    W(buf.str());

    // Hardness estimates
    sep();
    W("  Hardness estimates - Vickers Hv (GPa)\n");
    sep();
    buf.str(""); buf<<std::fixed<<std::setprecision(2);
    buf<<"  Chen et al. 2011  Hv = 2(k^2G)^0.585 - 3, k=G_H/K_H  :"
       <<std::setw(8)<<res.hardnessChen<<" GPa"
       <<"  [Intermetallics 19 (2011) 1275]\n"
       <<"  Tian et al. 2012  Hv = 0.92k^1.137 G_H^0.708       :"
       <<std::setw(8)<<res.hardnessTian<<" GPa"
       <<"  [IJRMHM 33 (2012) 93]\n"
       <<"  Teter     1998    Hv = 0.151 G_V                    :"
       <<std::setw(8)<<res.hardnessTeter<<" GPa"
       <<"  [MRS Bull. 23 (1998) 22]\n"
       <<"  Niu et al. 2012   Hv = G_H(1-2v)/3                 :"
       <<std::setw(8)<<res.hardnessNiu<<" GPa"
       <<"  [J. Phys.: CM 24 (2012) 405401]\n";
    W(buf.str());

    // Directional elastic properties
    sep();
    W("  Directional elastic properties\n");
    W("  Ref: Nye, Physical Properties of Crystals (1957)\n");
    sep();
    buf.str(""); buf<<std::fixed<<std::setprecision(3);
    buf<<"  Young's modulus E = 1/S_ii along principal axes:\n"
       <<"    E[100] = "<<std::setw(9)<<res.E_x<<" GPa"
       <<"    E[010] = "<<std::setw(9)<<res.E_y<<" GPa"
       <<"    E[001] = "<<std::setw(9)<<res.E_z<<" GPa\n";
    if (crystalFamily == "cubic") {
        buf<<"  Off-axis Young's moduli (cubic symmetry):\n"
           <<"    E[110] = "<<std::setw(9)<<res.E110<<" GPa"
           <<"    E[111] = "<<std::setw(9)<<res.E111<<" GPa\n";
    }
    buf<<"  Linear compressibility beta = sum_j S_ij (10^-3/GPa, dL/L per unit pressure):\n"
       <<std::setprecision(4)
       <<"    beta_x = "<<std::setw(9)<<res.beta_x*1e3
       <<"    beta_y = "<<std::setw(9)<<res.beta_y*1e3
       <<"    beta_z = "<<std::setw(9)<<res.beta_z*1e3<<"  (x10^-3 GPa^-1)\n";
    W(buf.str());

    // Sound velocities & thermal
    sep();
    W("  Sound velocities & thermal properties\n");
    sep();
    buf.str(""); buf<<std::fixed<<std::setprecision(1);
    buf<<"  v_L (longitudinal)    = "<<std::setw(10)<<res.vL<<" m/s\n"
       <<"  v_T (transverse)      = "<<std::setw(10)<<res.vT<<" m/s\n"
       <<"  v_m (mean/Debye)      = "<<std::setw(10)<<res.vM<<" m/s\n";
    buf<<"  Debye temp. Theta_D   = "<<std::setw(10)<<res.thetaDebye<<" K"
       <<"  [Anderson, J. Phys. Chem. Solids 24 (1963) 909]\n";
    buf<<std::setprecision(3)
       <<"  Gruneisen param. gamma = "<<std::setw(10)<<res.gruneisen
       <<"  [Belomestnykh & Tesleva, Tech. Phys. 49 (2004) 1304]\n";
    buf<<std::setprecision(4)
       <<"  kappa_min (CWP)       = "<<std::setw(10)<<res.kMinClarke<<" W/(m*K)"
       <<"  [Cahill, Watson & Pohl, PRB 46 (1992) 6131]\n"
       <<"  kappa_Slack (300K)    = "<<std::setw(10)<<res.kSlack   <<" W/(m*K)"
       <<"  [Slack, J. Phys. Chem. Solids 34 (1973) 321]\n";
    buf<<std::setprecision(1)
       <<"  Empirical T_m         = "<<std::setw(10)<<res.meltingTemp<<" K"
       <<"  [Fine, Brown & Marcus, Scr. Metall. 18 (1984) 951]\n";
    W(buf.str());

    // Mechanical stability
    sep();
    W("  Mechanical stability\n");
    W("  Ref: Born & Huang, Dynamical Theory of Crystal Lattices (1954);\n"
      "       Mouhat & Coudert, Phys. Rev. B 90 (2014) 224104\n");
    sep();
    if (res.mechanicallyStable) {
        W("  STABLE: all Born criteria satisfied and all C eigenvalues > 0\n");
    } else {
        W("  NOT STABLE:\n");
        for (const auto& n : res.stabilityNotes) W("    * "+n+"\n");
    }
    sep();
    W("\n");
}

}  // namespace qe
