#pragma once

#include <Eigen/Dense>

#include <array>
#include <string>
#include <vector>

namespace qe {

struct Atom {
    std::string symbol;
    Eigen::Vector3d fracPosition;
};

struct CifStructure {
    Eigen::Matrix3d cellAngstrom;
    std::vector<Atom> atoms;
    std::string spaceGroupName;
    int spaceGroupNumber = -1;
    std::vector<std::string> symOps;
};

struct BandData {
    std::vector<std::vector<double>> kByBand;
    std::vector<std::vector<double>> eByBand;
    // High-symmetry point markers: (k-distance, label). Label may be empty.
    std::vector<std::pair<double, std::string>> kLabelMarks;
};

// Projected band-structure data parsed from atomic_proj.xml (projwfc.x output)
struct AtomicProj {
    struct WfcInfo {
        int         atomnum = 0;  // 1-based atom index
        std::string elem;         // element name e.g. "Si"
        int         l = 0;        // angular momentum quantum number
    };
    int nk   = 0;
    int nbnd = 0;
    int nwfc = 0;
    std::vector<WfcInfo> wfcInfo;  // size = nwfc
    // weights[iwfc][ik][ibnd]  =  |<psi_nk | phi_iwfc>|^2
    std::vector<std::vector<std::vector<double>>> weights;
};

struct KPathNode {
    std::string label;
    Eigen::Vector3d k;
};

struct SymmetryKPath {
    std::string family;
    std::vector<KPathNode> nodes;
};

// ─── Elastic constants ────────────────────────────────────────────────────────

// Voigt notation index pairs: (1,1)=0 (2,2)=1 (3,3)=2 (2,3)=3 (1,3)=4 (1,2)=5
// Strain types used in the energy-strain method (one per independent C row).
struct StrainPattern {
    std::string name;           // e.g. "e1", "e4"
    std::array<double, 6> eta;  // Voigt strain vector (multiplied by δ at generation)
};

struct ElasticSetup {
    std::string crystalFamily;  // cubic, hexagonal, tetragonal, orthorhombic, …
    std::vector<StrainPattern> patterns;
    std::vector<double> deltas;  // strain magnitudes to sample (e.g. -0.04…+0.04)
    int nDeltas = 7;
};

struct ElasticResults {
    std::string crystalFamily;          // spglib-detected family
    Eigen::Matrix<double, 6, 6> C;  // stiffness matrix (GPa)
    Eigen::Matrix<double, 6, 6> S;  // compliance matrix (GPa⁻¹)

    // VRH averages
    double KV = 0, KR = 0, KH = 0;  // bulk modulus (GPa)
    double GV = 0, GR = 0, GH = 0;  // shear modulus (GPa)
    double EH = 0;                    // Young's modulus (GPa)
    double nuH = 0;                   // Poisson's ratio
    double AZ = 0;                    // Zener anisotropy (cubic)

    // Lame constants (Hill, GPa)
    double lambdaH = 0;  // first Lame parameter λ = K - 2G/3
    double muH = 0;      // second Lame parameter μ = G

    // Mechanical character
    double pughRatio = 0;      // K_H/G_H  (>1.75 ductile, <1.75 brittle)
    double cauchyPressure = 0; // C12 - C44 (GPa)

    // Anisotropy indices
    double AU = 0;  // universal anisotropy (Ranganathan & Ostoja-Starzewski, 2008)
    double AB = 0;  // percent bulk anisotropy  (%)
    double AG = 0;  // percent shear anisotropy (%)

    // Hardness estimates (GPa)
    double hardnessChen  = 0;  // Chen et al., Intermetallics 19 (2011) 1275
    double hardnessTian  = 0;  // Tian et al., IJRMHM 33 (2012) 93
    double hardnessTeter = 0;  // Teter, MRS Bull. 23 (1998) 22
    double hardnessNiu   = 0;  // Niu et al., J. Phys.: CM 24 (2012) 405401

    // Sound velocities (m/s)
    double vL = 0, vT = 0, vM = 0;  // longitudinal, transverse, mean
    double vDebye = 0;               // Debye velocity

    // Debye temperature (K)
    double thetaDebye = 0;

    // Thermal / acoustic properties
    double gruneisen   = 0;  // acoustic Grüneisen (Poisson's ratio route)
    double kMinClarke  = 0;  // min. thermal cond., Cahill-Watson-Pohl (W/m/K)
    double kSlack      = 0;  // Slack lattice thermal cond. at 300K (W/m/K)
    double meltingTemp = 0;  // empirical Tm, Fine-Brown-Marcus (K)

    // Directional Young's moduli (GPa) from compliance
    double E_x = 0, E_y = 0, E_z = 0;  // along [100],[010],[001]  (all systems)
    double E110 = 0, E111 = 0;          // along [110],[111]        (cubic only)

    // Linear compressibility (1/GPa) under hydrostatic pressure
    double beta_x = 0, beta_y = 0, beta_z = 0;

    // Per-plane shear anisotropy (Chung & Buessem, J. Appl. Phys. 38 (1967) 2535)
    // A = 1 for isotropic; {100},{010},{001} planes respectively
    double A1_plane = 0, A2_plane = 0, A3_plane = 0;

    // Stability
    bool mechanicallyStable = false;
    std::vector<std::string> stabilityNotes;
};

}  // namespace qe
