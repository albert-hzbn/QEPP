// QHA elastic constants: I/O layer.
// Pre-processing: generate volume-scaled inputs + elastic strain inputs + DFPT inputs.
// Post-processing: read summary, run elastic post-processing per volume,
//                  read phonon DOS, compute temperature-dependent C_ij(T).

#include "qe/qha_elastic.hpp"
#include "qe/elastic.hpp"
#include "qe/qha.hpp"
#include "qe/dfpt.hpp"
#include "qe/struct.hpp"
#include "qe/utils.hpp"

#include <Eigen/Dense>

#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <numeric>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace qe {

// ── Internal helpers ──────────────────────────────────────────────────────────

static std::string qhe_vol_dir_name(int idx, int total) {
    const int w = (total < 100) ? 2 : 3;
    std::ostringstream oss;
    oss << "v" << std::setw(w) << std::setfill('0') << (idx + 1);
    return oss.str();
}

// Resolve a path that may be relative to a base directory.
static std::string resolve_path(const std::string& base, const std::string& path) {
    if (path.empty()) return path;
    if (path[0] == '/') return path;  // already absolute
    return (fs::path(base) / path).string();
}

static std::string shell_quote(const std::string& text) {
    std::string out = "'";
    for (char c : text) {
        if (c == '\'') out += "'\\''";
        else out.push_back(c);
    }
    out += "'";
    return out;
}

static int run_shell_command(const std::string& cmd) {
    const int rc = std::system(cmd.c_str());
    if (rc == -1) return rc;
    return rc;
}

static bool contains_job_done(const fs::path& path) {
    std::ifstream in(path);
    if (!in.is_open()) return false;
    std::string line;
    while (std::getline(in, line))
        if (line.find("JOB DONE") != std::string::npos) return true;
    return false;
}

static double read_total_energy_from_output(const fs::path& outPath) {
    std::ifstream in(outPath);
    if (!in.is_open())
        throw std::runtime_error("Could not open SCF output: " + outPath.string());

    double energy = 0.0;
    bool found = false;
    std::string line;
    while (std::getline(in, line)) {
        if (line.find('!') == std::string::npos || line.find("total energy") == std::string::npos)
            continue;
        std::istringstream ss(line);
        std::string tmp;
        ss >> tmp >> tmp >> tmp >> tmp >> energy;
        found = true;
    }
    if (!found)
        throw std::runtime_error("No total energy found in: " + outPath.string());
    return energy;
}

static double resolve_static_energy(const std::string& energyToken,
                                    const fs::path& scfAbs) {
    double energy = 0.0;
    if (try_parse_double(energyToken, energy))
        return energy;

    const fs::path volumeDir = scfAbs.parent_path();
    const fs::path qeOut = volumeDir / "qe.out";
    if (fs::exists(qeOut)) return read_total_energy_from_output(qeOut);

    const fs::path stemOut = volumeDir / (scfAbs.stem().string() + ".out");
    if (fs::exists(stemOut)) return read_total_energy_from_output(stemOut);

    throw std::runtime_error("Energy column is not numeric and no SCF output was found for " + scfAbs.string());
}

// ── Pre-processing ─────────────────────────────────────────────────────────────

void qha_elastic_generate_inputs(const std::string& qeInputPath,
                                   int    nVolumes,
                                   double rangePercent,
                                   const std::string& outDirArg,
                                   int    nDeltas,
                                   double maxDelta,
                                   const DfptOptions& dfptOpts) {
    if (nVolumes < 4)
        throw std::runtime_error("QHA elastic requires at least 4 volume points.");
    if (rangePercent <= 0.0 || rangePercent >= 100.0)
        throw std::runtime_error("rangePercent must be in (0, 100).");
    if (nDeltas < 3 || nDeltas % 2 == 0)
        throw std::runtime_error("nDeltas must be odd and >= 3.");
    if (maxDelta <= 0.0 || maxDelta >= 0.2)
        throw std::runtime_error("maxDelta must be in (0, 0.2).");

    // Parse equilibrium structure
    const StructInfo info = parse_struct_from_qe_input(qeInputPath);
    if (info.atoms.empty() || std::abs(info.cellAngst.determinant()) < 1e-12)
        throw std::runtime_error("Could not parse structure from: " + qeInputPath);

    const double V0    = std::abs(info.cellAngst.determinant());
    const int    natom = info.nAtoms;
    const double half  = 0.5 * rangePercent / 100.0;
    const double dFrac = (nVolumes > 1) ?
        (2.0 * half / static_cast<double>(nVolumes - 1)) : 0.0;

    const std::string outDir = outDirArg.empty() ?
        (stem_from_path(qeInputPath) + "_qha_elastic") : outDirArg;
    fs::create_directories(outDir);

    const auto origLines = load_lines(qeInputPath);
    std::string qePrefix = extract_quoted_assignment(origLines, "prefix");
    if (qePrefix.empty()) qePrefix = stem_from_path(qeInputPath);
    const std::string scfStem = stem_from_path(qeInputPath);

    // Summary file
    const std::string summaryPath = outDir + "/qha_elastic_summary.in";
    std::ofstream fsum(summaryPath);
    if (!fsum.is_open())
        throw std::runtime_error("Could not create: " + summaryPath);

    fsum << "# QHA elastic summary file generated by qepp qha_elastic -pre\n";
    fsum << "#\n";
    fsum << "# Format: volume(Ang^3)   energy(Ry)   scf_template   elastic_dir   phonon_dos_path\n";
    fsum << "#\n";
    fsum << "# Workflow for each v*/ directory:\n";
    fsum << "#   1. pw.x  < " << scfStem << ".in                    (SCF at this volume)\n";
    fsum << "#   2. for d in elastic/*/*/; do (cd $d && pw.x < "
         << qePrefix << ".in > " << qePrefix << ".out); done\n";
    fsum << "#   3. ph.x  < dfpt/ph.in  > dfpt/ph.out\n";
    fsum << "#   4. q2r.x < dfpt/q2r.in\n";
    fsum << "#   5. matdyn.x < dfpt/matdyn_dos.in\n";
    fsum << "#      (matdyn.x writes " << qePrefix << ".phonon.dos to its CWD)\n";
    fsum << "#\n";
    fsum << "# Fill in the energy(Ry) column from the SCF output of each v*/:\n";
    fsum << "#   grep '!.*total energy' v01/" << scfStem << ".out | tail -1 | awk '{print $5}'\n";
    fsum << "#\n";
    fsum << "# Then run:\n";
    fsum << "#   qepp qha_elastic -post " << summaryPath
         << " [--tmin T] [--tmax T] [--dt T]\n";
    fsum << "#\n";
    fsum << std::left
         << std::setw(18) << "# vol(Ang^3)"
         << std::setw(24) << "energy(Ry)"
         << std::setw(24) << "scf_template"
         << std::setw(24) << "elastic_dir"
         << "phonon_dos_path\n";

    int totalElasticInputs = 0;
    int firstElasticPatterns = 0;  // printed once for the workflow message

    for (int i = 0; i < nVolumes; ++i) {
        const double scale   = 1.0 - half + i * dFrac;
        const double s       = std::cbrt(scale);          // linear scale
        const Eigen::Matrix3d newCell = info.cellAngst * s;
        const double Vi = std::abs(newCell.determinant());

        const std::string vdir  = outDir + "/" + qhe_vol_dir_name(i, nVolumes);
        const std::string scfIn = vdir + "/" + scfStem + ".in";
        fs::create_directories(vdir);

        // Write volume-scaled SCF input (uses same helper logic as qha.cpp)
        {
            // We call qha_generate_volumes logic inlined here so we can also
            // create elastic directories before generating their inputs.
            // Reuse the public utils: load_lines + extract_quoted_assignment.
            std::ostringstream cellBlock;
            cellBlock << std::fixed << std::setprecision(10);
            bool inSystem = false, inCell = false;
            int  cellRow  = 0;
            bool wroteCell = false;
            std::ofstream fout(scfIn);
            if (!fout.is_open())
                throw std::runtime_error("Cannot create: " + scfIn);

            for (const auto& raw : origLines) {
                const std::string t  = trim(raw);
                const std::string lo = to_lower(t);

                if (lo.rfind("&system", 0) == 0) inSystem = true;
                if (inSystem && t == "/")         inSystem = false;

                if (inSystem) {
                    if (lo.find("ibrav") != std::string::npos && lo.find('=') != std::string::npos) {
                        fout << "  ibrav = 0,\n"; continue;
                    }
                    if (lo.find("celldm") != std::string::npos && lo.find('=') != std::string::npos) {
                        continue;
                    }
                }
                if (lo.find("outdir") != std::string::npos && lo.find('=') != std::string::npos) {
                    fout << "  outdir = './tmp'\n"; continue;
                }
                if (lo.rfind("cell_parameters", 0) == 0) {
                    inCell = true; cellRow = 0;
                    fout << "CELL_PARAMETERS angstrom\n";
                    wroteCell = true; continue;
                }
                if (inCell && cellRow < 3) {
                    std::istringstream ss(t); double x, y, z;
                    if (ss >> x >> y >> z) {
                        fout << "  " << std::setw(18) << std::fixed << std::setprecision(10)
                             << newCell(cellRow, 0) << "  " << newCell(cellRow, 1)
                             << "  " << newCell(cellRow, 2) << "\n";
                        if (++cellRow == 3) inCell = false;
                        continue;
                    }
                }
                fout << raw << "\n";
            }
            if (!wroteCell) {
                fout << "CELL_PARAMETERS angstrom\n";
                for (int r = 0; r < 3; ++r)
                    fout << "  " << std::setw(18) << std::fixed << std::setprecision(10)
                         << newCell(r,0) << "  " << newCell(r,1)
                         << "  " << newCell(r,2) << "\n";
            }
        }

        // Elastic strain inputs for this volume
        const std::string elasticDir = vdir + "/elastic";
        try {
            generate_elastic_inputs(scfIn, elasticDir, nDeltas, maxDelta);
        } catch (const std::exception& e) {
            std::cerr << "  WARNING: elastic inputs for " << vdir
                      << " failed: " << e.what() << "\n";
        }

        // DFPT phonon inputs for this volume
        const std::string dfptDir = vdir + "/dfpt";
        try {
            generate_phonon_inputs(scfIn, dfptDir, dfptOpts, /*verbose=*/false);
        } catch (const std::exception& e) {
            std::cerr << "  WARNING: phonon inputs for " << vdir
                      << " failed: " << e.what() << "\n";
        }

        // phonon.dos is written by matdyn.x to the volume directory
        const std::string volDirRel = qhe_vol_dir_name(i, nVolumes);
        fsum << std::fixed << std::setprecision(6) << std::left
             << std::setw(18) << Vi
             << std::setw(24) << "TODO_fill_energy"
             << std::setw(24) << (volDirRel + "/" + scfStem + ".in")
             << std::setw(24) << (volDirRel + "/elastic")
             << volDirRel << "/" << qePrefix << ".phonon.dos\n";

        if (i == 0) {
            // Count elastic inputs from the first volume for summary message
            try {
                std::ifstream mf(elasticDir + "/elastic_setup.dat");
                int nd = 0; std::string key; std::string line;
                while (std::getline(mf, line)) {
                    std::istringstream ss(line); ss >> key;
                    if (key == "ndeltas") ss >> nd;
                }
                // Rough pattern count from directory listing
                int npat = 0;
                for (const auto& e : fs::directory_iterator(elasticDir))
                    if (e.is_directory()) ++npat;
                firstElasticPatterns = npat;
                totalElasticInputs += npat * nd;
            } catch (...) {}
        } else {
            totalElasticInputs += firstElasticPatterns * nDeltas;
        }

        std::cout << "  Created " << vdir << "  (V = " << std::fixed
                  << std::setprecision(4) << Vi << " Å³, scale = "
                  << std::setprecision(4) << scale << ")\n";
    }

    (void)natom;

    std::cout << "\nOutput directory : " << outDir << "\n";
    std::cout << "Summary file     : " << summaryPath << "\n";
    std::cout << "\nTotal calculations needed per volume:\n";
    std::cout << "  1 × SCF (volume energy)\n";
    std::cout << "  " << firstElasticPatterns * nDeltas
              << " × SCF (elastic strain deformations)\n";
    std::cout << "  1 × ph.x + q2r.x + matdyn.x (phonon DOS)\n";
    std::cout << "\nTotal across all " << nVolumes << " volumes: "
              << nVolumes * (1 + firstElasticPatterns * nDeltas + 1)
              << " SCF + " << nVolumes << " phonon runs\n";
    std::cout << "\nFill in energy(Ry) for each volume, then run:\n";
    std::cout << "  qepp qha_elastic -post " << summaryPath
              << " --tmin 0 --tmax 1500 --dt 10\n";
}

void qha_elastic_run_dataset(const std::string& datasetDir,
                             const QeParallelOptions& parallel,
                             const std::set<std::string>& excludeVolumes) {
    if (parallel.np < 1)
        throw std::runtime_error("--np must be >= 1.");

    const fs::path dataset = fs::path(datasetDir);
    if (!fs::is_directory(dataset))
        throw std::runtime_error("Dataset directory not found: " + datasetDir);

    std::vector<fs::path> volumes;
    for (const auto& entry : fs::directory_iterator(dataset)) {
        if (!entry.is_directory()) continue;
        const std::string name = entry.path().filename().string();
        if (name.size() >= 2 && name[0] == 'v' && std::isdigit(static_cast<unsigned char>(name[1])))
            volumes.push_back(entry.path());
    }
    std::sort(volumes.begin(), volumes.end());
    if (volumes.empty())
        throw std::runtime_error("No v*/ directories found under: " + datasetDir);

    const std::string launch = qe_parallel_args(parallel);

    for (const auto& vdir : volumes) {
        const std::string vname = vdir.filename().string();
        if (excludeVolumes.count(vname)) {
            std::cout << "===== " << vname << " =====\n";
            std::cout << "  skipped by --exclude\n";
            continue;
        }

        std::cout << "===== " << vname << " =====\n";
        fs::path scfIn;
        for (const auto& entry : fs::directory_iterator(vdir)) {
            if (entry.is_regular_file() && entry.path().extension() == ".in") {
                scfIn = entry.path();
                break;
            }
        }
        if (scfIn.empty())
            throw std::runtime_error("No top-level SCF input found in: " + vdir.string());

        const fs::path qeOut = vdir / "qe.out";
        if (!contains_job_done(qeOut)) {
            std::cout << "  SCF: running...\n";
            const std::string cmd = "cd " + shell_quote(vdir.string()) +
                " && " + launch + " pw.x -input " + shell_quote(scfIn.filename().string()) +
                " > qe.out 2>&1";
            run_shell_command(cmd);
            if (!contains_job_done(qeOut)) {
                std::cout << "  SCF FAILED\n";
                continue;
            }
            std::cout << "  SCF: done\n";
        } else {
            std::cout << "  SCF: done\n";
        }

        const fs::path elasticDir = vdir / "elastic";
        for (const auto& patt : fs::directory_iterator(elasticDir)) {
            if (!patt.is_directory()) continue;
            for (const auto& strain : fs::directory_iterator(patt.path())) {
                if (!strain.is_directory()) continue;
                fs::path infile;
                for (const auto& f : fs::directory_iterator(strain.path())) {
                    if (f.is_regular_file() && f.path().extension() == ".in") {
                        infile = f.path();
                        break;
                    }
                }
                if (infile.empty()) continue;
                const fs::path outfile = infile.parent_path() / (infile.stem().string() + ".out");
                if (contains_job_done(outfile)) continue;
                const std::string cmd = "cd " + shell_quote(strain.path().string()) +
                    " && " + launch + " pw.x < " + shell_quote(infile.filename().string()) +
                    " > " + shell_quote(outfile.filename().string()) + " 2>&1";
                run_shell_command(cmd);
            }
        }

        const fs::path tmpPh = vdir / "tmp" / "_ph0";
        fs::create_directories(tmpPh);
        const fs::path phOut = vdir / "dfpt" / "ph.out";
        if (!contains_job_done(phOut)) {
            std::cout << "  DFPT: running...\n";
            const std::string cmd = "cd " + shell_quote(vdir.string()) +
                " && " + launch + " ph.x -input dfpt/ph.in > dfpt/ph.out 2>&1";
            run_shell_command(cmd);
            if (!contains_job_done(phOut)) {
                std::cout << "  DFPT FAILED\n";
                continue;
            }
            std::cout << "  DFPT: done\n";
        } else {
            std::cout << "  DFPT: done\n";
        }

        const fs::path q2rOut = vdir / "dfpt" / "q2r.out";
        if (!contains_job_done(q2rOut)) {
            std::cout << "  q2r: running...\n";
            const std::string cmd = "cd " + shell_quote(vdir.string()) +
                " && " + launch + " q2r.x < dfpt/q2r.in > dfpt/q2r.out 2>&1";
            run_shell_command(cmd);
            if (!contains_job_done(q2rOut)) {
                std::cout << "  q2r FAILED\n";
                continue;
            }
            std::cout << "  q2r: done\n";
        } else {
            std::cout << "  q2r: done\n";
        }

        const fs::path matdynOut = vdir / "dfpt" / "matdyn_dos.out";
        if (!contains_job_done(matdynOut)) {
            std::cout << "  matdyn: running...\n";
            const std::string cmd = "cd " + shell_quote(vdir.string()) +
                " && " + launch + " matdyn.x < dfpt/matdyn_dos.in > dfpt/matdyn_dos.out 2>&1";
            run_shell_command(cmd);
            if (!contains_job_done(matdynOut)) {
                std::cout << "  matdyn FAILED\n";
                continue;
            }
            std::cout << "  matdyn: done\n";
        } else {
            std::cout << "  matdyn: done\n";
        }

        std::cout << "  " << vname << ": COMPLETE\n";
    }
}

// ── Post-processing ────────────────────────────────────────────────────────────

QhaElasticResult read_and_compute_qha_elastic(const std::string& summaryPath,
                                               const std::vector<double>& tempsIn,
                                               const std::set<std::string>& excludeVolumes) {
    const fs::path baseDir = fs::path(summaryPath).parent_path();

    std::ifstream fin(summaryPath);
    if (!fin.is_open())
        throw std::runtime_error("Cannot open summary file: " + summaryPath);

    std::vector<QhaElasticVolumePoint> elasticVols;
    std::vector<QhaVolumePoint>        phonVols;

    std::string line;
    while (std::getline(fin, line)) {
        const std::string t = trim(line);
        if (t.empty() || t[0] == '#') continue;

        std::istringstream ss(t);
        double vol;
        std::string energyToken;
        std::string scfTemplate, elasticDir, dosPath;
        if (!(ss >> vol >> energyToken >> scfTemplate >> elasticDir >> dosPath))
            throw std::runtime_error("Malformed line in summary file: " + line);

        const std::string scfAbs     = resolve_path(baseDir.string(), scfTemplate);
        const std::string elasticAbs = resolve_path(baseDir.string(), elasticDir);
        const std::string dosAbs     = resolve_path(baseDir.string(), dosPath);
        const std::string volumeName = fs::path(scfTemplate).parent_path().filename().string();
        if (!volumeName.empty() && excludeVolumes.count(volumeName))
            continue;
        const double energy = resolve_static_energy(energyToken, fs::path(scfAbs));

        // --- Elastic C_ij at this volume ---
        std::cout << "  Computing elastic constants for V = " << std::fixed
                  << std::setprecision(4) << vol << " Å³ ...\n";
        ElasticResults er;
        try {
            er = compute_elastic_properties(scfAbs, elasticAbs);
        } catch (const std::exception& e) {
            throw std::runtime_error("Elastic post-processing failed for volume " +
                std::to_string(vol) + " Å³: " + e.what());
        }

        QhaElasticVolumePoint ep;
        ep.volumeAng3 = vol;
        ep.energyRy   = energy;
        ep.elastic    = er;
        elasticVols.push_back(ep);

        // --- Phonon free energy at this volume ---
        std::cout << "  Reading phonon DOS for V = " << std::fixed
                  << std::setprecision(4) << vol << " Å³ ...\n";
        QhaVolumePoint qp;
        try {
            qp = compute_qha_volume_from_dos(dosAbs, vol, energy, tempsIn);
        } catch (const std::exception& e) {
            throw std::runtime_error("Phonon DOS processing failed for volume " +
                std::to_string(vol) + " Å³: " + e.what());
        }
        phonVols.push_back(qp);
    }

    if (elasticVols.size() < 4)
        throw std::runtime_error(
            "Need at least 4 volume entries in summary file, found " +
            std::to_string(elasticVols.size()));

    // Sort both arrays by volume (ascending)
    std::vector<size_t> idx(elasticVols.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(),
              [&](size_t a, size_t b){ return elasticVols[a].volumeAng3 < elasticVols[b].volumeAng3; });

    std::vector<QhaElasticVolumePoint> eVols(elasticVols.size());
    std::vector<QhaVolumePoint>        pVols(phonVols.size());
    for (size_t k = 0; k < idx.size(); ++k) {
        eVols[k] = elasticVols[idx[k]];
        pVols[k] = phonVols[idx[k]];
    }

    std::cout << "\nRunning QHA elastic computation ...\n";
    return compute_qha_elastic(eVols, pVols);
}

// ── Report ─────────────────────────────────────────────────────────────────────

void write_qha_elastic_report(const QhaElasticResult& result,
                               const std::string& outPrefix) {
    const std::string fname = outPrefix + ".qha_elastic.txt";

    auto write_report = [&](std::ostream& os) {
        os << "# Temperature-Dependent Elastic Constants (QHA)\n";
        os << "# Crystal family: " << result.crystalFamily << "\n";
        os << "# N volumes used: " << result.nvolumes << "\n";
        os << "#\n";
        os << "# Method: Quasi-harmonic approximation\n";
        os << "#   C_ij(T) = C_ij^{QS}(V(T)) + ΔC_ij^{ph}(T)\n";
        os << "#   Quasi-static: C_ij polynomial interpolated to QHA volume V(T)\n";
        os << "#   Phonon correction: B^{ph}=V·d²F_vib/dV² added to normal-normal entries\n";
        os << "#\n";

        // Column header: T, V, alpha, then all 21 unique C_ij elements,
        //                then bulk/shear/young/poisson from VRH, then phonon B
        os << "# " << std::left
           << std::setw(10) << "T(K)"
           << std::setw(14) << "V(Ang^3)"
           << std::setw(14) << "alpha(1/K)";

        // Print unique upper-triangle C_ij labels
        for (int r = 0; r < 6; ++r)
            for (int c = r; c < 6; ++c) {
                std::string lbl = "C" + std::to_string(r+1) + std::to_string(c+1) + "(GPa)";
                os << std::setw(14) << lbl;
            }
        os << std::setw(12) << "KH(GPa)"
           << std::setw(12) << "GH(GPa)"
           << std::setw(12) << "EH(GPa)"
           << std::setw(10) << "nuH"
           << std::setw(14) << "B_QS(GPa)"
           << std::setw(14) << "B_ph(GPa)"
           << "\n";

        os << std::fixed;
        for (const auto& pt : result.thermal) {
            os << std::setw(12) << std::setprecision(1)  << pt.temperature
               << std::setw(14) << std::setprecision(6)  << pt.volumeAng3
               << std::setw(14) << std::setprecision(6)  << pt.alpha;

            for (int r = 0; r < 6; ++r)
                for (int c = r; c < 6; ++c)
                    os << std::setw(14) << std::setprecision(4) << pt.C_total(r,c);

            os << std::setw(12) << std::setprecision(4) << pt.KH
               << std::setw(12) << std::setprecision(4) << pt.GH
               << std::setw(12) << std::setprecision(4) << pt.EH
               << std::setw(10) << std::setprecision(4) << pt.nuH
               << std::setw(14) << std::setprecision(4) << pt.bulkModStatic
               << std::setw(14) << std::setprecision(4) << pt.bulkModPhonon
               << "\n";
        }
    };

    // stdout
    std::cout << "\n";
    write_report(std::cout);

    // file
    std::ofstream fout(fname);
    if (fout.is_open()) {
        write_report(fout);
        std::cout << "\nData written to: " << fname << "\n";
    } else {
        std::cerr << "Warning: could not write report to " << fname << "\n";
    }

    // Also print a concise summary table to stdout
    std::cout << "\n"
              << std::string(80, '=') << "\n"
              << "  T-dependent elastic constants summary"
              << " (" << result.crystalFamily << ")\n"
              << std::string(80, '-') << "\n"
              << std::left
              << std::setw(8)  << "T(K)"
              << std::setw(12) << "V(Ang^3)"
              << std::setw(10) << "alpha(1/K)"
              << std::setw(10) << "KH(GPa)"
              << std::setw(10) << "GH(GPa)"
              << std::setw(10) << "EH(GPa)"
              << std::setw(8)  << "nuH"
              << std::setw(10) << "B_ph\n"
              << std::string(80, '-') << "\n"
              << std::fixed;

    for (const auto& pt : result.thermal) {
        std::cout << std::setw(8)  << std::setprecision(0) << pt.temperature
                  << std::setw(12) << std::setprecision(4) << pt.volumeAng3
                  << std::setw(10) << std::setprecision(2) << pt.alpha
                  << std::setw(10) << std::setprecision(2) << pt.KH
                  << std::setw(10) << std::setprecision(2) << pt.GH
                  << std::setw(10) << std::setprecision(2) << pt.EH
                  << std::setw(8)  << std::setprecision(4) << pt.nuH
                  << std::setw(10) << std::setprecision(2) << pt.bulkModPhonon
                  << "\n";
    }
    std::cout << std::string(80, '=') << "\n";
}

}  // namespace qe
