#include "qe/cli/commands.hpp"

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "qe/dfpt.hpp"
#include "qe/help.hpp"
#include "qe/phonon.hpp"
#include "qe/utils.hpp"

namespace qe {

// ─── Shared helpers ───────────────────────────────────────────────────────────

// Parse comma-separated labels from a string into a vector.
static std::vector<std::string> parse_labels(const std::string& lst) {
    std::vector<std::string> out;
    std::string tok;
    for (char c : lst) {
        if (c == ',') { if (!tok.empty()) out.push_back(tok); tok.clear(); }
        else tok += c;
    }
    if (!tok.empty()) out.push_back(tok);
    return out;
}

// Build a temperature grid from tMin to tMax (inclusive) with step dt.
static std::vector<double> make_temp_grid(double tMin, double tMax, double dt) {
    std::vector<double> temps;
    for (double T = tMin; T <= tMax + 1e-6; T += dt)
        temps.push_back(T);
    return temps;
}

// ── qepp phonon -pre <scf.in>
//       [--nq NQ1 NQ2 NQ3] [--nq-dos NQ1 NQ2 NQ3] [--epsil] [outdir] ─────────
int handle_phonon_pre_mode(int argc, char** argv, int s) {
    if (argc < 3 + s) {
        print_help_command(argv[0], "phonon", "-pre");
        return 1;
    }

    const std::string scfPath = argv[2 + s];
    std::string outDir;
    bool outDirSet = false;
    DfptOptions opts;

    for (int i = 3 + s; i < argc; ++i) {
        const std::string arg = argv[i];
        auto next3Int = [&](int& a, int& b, int& c2, const char* flag) {
            if (i + 3 >= argc)
                throw std::runtime_error(std::string(flag) + " requires 3 integers.");
            a = std::stoi(argv[++i]);
            b = std::stoi(argv[++i]);
            c2 = std::stoi(argv[++i]);
        };
        if (arg == "--nq") {
            next3Int(opts.nq1, opts.nq2, opts.nq3, "--nq");
        } else if (arg == "--nq-dos") {
            next3Int(opts.nqDos1, opts.nqDos2, opts.nqDos3, "--nq-dos");
        } else if (arg == "--epsil") {
            opts.epsil = true;
        } else if (arg == "--asr") {
            if (i + 1 >= argc) throw std::runtime_error("--asr requires a value.");
            opts.asr = argv[++i];
        } else if (arg == "--outdir") {
            if (i + 1 >= argc) throw std::runtime_error("--outdir requires a directory path.");
            outDir = argv[++i];
            outDirSet = true;
        } else if (!outDirSet && arg.rfind("--", 0) != 0) {
            outDir = arg;
            outDirSet = true;
        } else {
            throw std::runtime_error("Unknown argument for 'phonon -pre': " + arg);
        }
    }

    if (outDir.empty())
        outDir = stem_from_path(scfPath) + "_phonon";
    generate_phonon_inputs(scfPath, outDir, opts);
    return 0;
}

// ── qepp phonon -dos <dos_file> [--pdos <file>] [--labels a,b] [outprefix] ────
// Accepts both phonopy total_dos.dat and matdyn.x DOS output.
// Format is auto-detected from the file header (cm-1 marker → matdyn).
int handle_phonon_dos_mode(int argc, char** argv, int s) {
    if (argc < 3 + s) {
        print_help_command(argv[0], "phonon", "-dos");
        return 1;
    }

    const std::string dosPath = argv[2 + s];
    std::string projDosPath;
    std::vector<std::string> projLabels;
    std::string outPrefix = stem_from_path(dosPath);
    bool outPrefixSet = false;

    for (int i = 3 + s; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--pdos") {
            if (i + 1 >= argc) throw std::runtime_error("--pdos requires a file path.");
            projDosPath = argv[++i];
        } else if (arg == "--labels") {
            if (i + 1 >= argc) throw std::runtime_error("--labels requires a comma-separated list.");
            projLabels = parse_labels(argv[++i]);
        } else if (!outPrefixSet && arg.rfind("--", 0) != 0) {
            outPrefix = arg;
            outPrefixSet = true;
        } else {
            throw std::runtime_error("Unknown argument for 'phonon -dos': " + arg);
        }
    }

    PhononDosData dosData;
    if (is_matdyn_dos_format(dosPath)) {
        dosData = parse_matdyn_dos(dosPath);
    } else {
        dosData = parse_phonon_dos(dosPath, projDosPath, projLabels);
    }
    write_phonon_dos_plot(dosData, outPrefix);
    return 0;
}

// ── qepp phonon -band <file> [--labels G,X,M,G] [--qlabels <file>] [outprefix]
// Accepts both phonopy band.yaml and matdyn.x &plot frequency output.
// Format is auto-detected: .yaml extension or "nqpoint:" header → phonopy;
// "&plot" header or .freq extension → matdyn.
int handle_phonon_band_mode(int argc, char** argv, int s) {
    if (argc < 3 + s) {
        print_help_command(argv[0], "phonon", "-band");
        return 1;
    }

    const std::string bandPath = argv[2 + s];
    std::string outPrefix = stem_from_path(bandPath);
    std::string qLabelsFile;
    std::vector<std::string> manualLabels;
    bool outPrefixSet = false;

    for (int i = 3 + s; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--labels") {
            if (i + 1 >= argc) throw std::runtime_error("--labels requires a comma-separated list.");
            manualLabels = parse_labels(argv[++i]);
        } else if (arg == "--qlabels") {
            if (i + 1 >= argc) throw std::runtime_error("--qlabels requires a file path.");
            qLabelsFile = argv[++i];
        } else if (!outPrefixSet && arg.rfind("--", 0) != 0) {
            outPrefix = arg;
            outPrefixSet = true;
        } else {
            throw std::runtime_error("Unknown argument for 'phonon -band': " + arg);
        }
    }

    PhononBandData bandData;
    if (is_matdyn_freq_format(bandPath)) {
        bandData = parse_matdyn_freq(bandPath, qLabelsFile);
    } else {
        bandData = parse_phonon_band_yaml(bandPath);
    }

    // Apply manual labels if provided: distribute evenly across labelMarks
    if (!manualLabels.empty()) {
        if (manualLabels.size() == bandData.labelMarks.size()) {
            for (size_t k = 0; k < manualLabels.size(); ++k)
                bandData.labelMarks[k].second = manualLabels[k];
        } else {
            // Redistribute: assume first and last, then interpolate positions
            bandData.labelMarks.clear();
            const double dTotal = bandData.distances.empty() ? 1.0 : bandData.distances.back();
            const int nSeg = static_cast<int>(manualLabels.size()) - 1;
            if (nSeg > 0) {
                for (size_t k = 0; k < manualLabels.size(); ++k) {
                    const double d = dTotal * static_cast<double>(k) / nSeg;
                    bandData.labelMarks.push_back({d, manualLabels[k]});
                }
            }
        }
    }

    write_phonon_band_plot(bandData, outPrefix);
    return 0;
}

// ── qepp phonon -ha <dos_file>
//       [--tmin T] [--tmax T] [--dt T] [--natom N] [outprefix] ─────────────────
// Accepts both phonopy total_dos.dat and matdyn.x DOS output (auto-detected).
int handle_phonon_ha_mode(int argc, char** argv, int s) {
    if (argc < 3 + s) {
        print_help_command(argv[0], "phonon", "-ha");
        return 1;
    }

    const std::string dosPath = argv[2 + s];
    std::string outPrefix = stem_from_path(dosPath);
    bool outPrefixSet = false;

    double tMin  = 50.0;
    double tMax  = 1000.0;
    double dt    = 50.0;
    int    natom = 1;

    for (int i = 3 + s; i < argc; ++i) {
        const std::string arg = argv[i];
        auto nextDouble = [&](const char* flag) -> double {
            if (i + 1 >= argc)
                throw std::runtime_error(std::string(flag) + " requires a value.");
            return std::stod(argv[++i]);
        };
        auto nextInt = [&](const char* flag) -> int {
            if (i + 1 >= argc)
                throw std::runtime_error(std::string(flag) + " requires a value.");
            return std::stoi(argv[++i]);
        };
        if      (arg == "--tmin")  tMin  = nextDouble("--tmin");
        else if (arg == "--tmax")  tMax  = nextDouble("--tmax");
        else if (arg == "--dt")    dt    = nextDouble("--dt");
        else if (arg == "--natom") natom = nextInt("--natom");
        else if (!outPrefixSet && arg.rfind("--", 0) != 0) {
            outPrefix = arg;
            outPrefixSet = true;
        } else {
            throw std::runtime_error("Unknown argument for 'phonon -ha': " + arg);
        }
    }

    if (tMax <= tMin || dt <= 0.0)
        throw std::runtime_error("Invalid temperature range: check --tmin/--tmax/--dt.");
    if (natom < 1)
        throw std::runtime_error("--natom must be >= 1.");

    PhononDosData dosData;
    if (is_matdyn_dos_format(dosPath)) {
        dosData = parse_matdyn_dos(dosPath);
    } else {
        dosData = parse_phonon_dos(dosPath);
    }

    const auto haResult = compute_harmonic_approximation(dosData, make_temp_grid(tMin, tMax, dt), natom);
    write_ha_report(haResult, outPrefix);
    write_ha_plot(haResult, outPrefix);
    return 0;
}

// ── qepp phonon -post <prefix>
//       [--dos] [--band] [--ha] [--natom N] [--tmin T] [--tmax T] [--dt T]
//       [--labels G,X,...] [--qlabels <file>] [outprefix] ─────────────────────
// Post-processes matdyn.x output: reads <prefix>.phonon.dos and/or
// <prefix>.phonon_band.freq (the default output names from phonon -pre).
// By default runs all three (DOS, band, HA); use flags to select individually.
int handle_phonon_post_mode(int argc, char** argv, int s) {
    if (argc < 3 + s) {
        print_help_command(argv[0], "phonon", "-post");
        return 1;
    }

    const std::string prefix = argv[2 + s];
    std::string outPrefix = prefix;
    bool outPrefixSet = false;

    bool doDos  = false;
    bool doBand = false;
    bool doHa   = false;

    double tMin  = 50.0;
    double tMax  = 1000.0;
    double dt    = 50.0;
    int    natom = 1;
    std::string qLabelsFile;
    std::vector<std::string> manualLabels;

    for (int i = 3 + s; i < argc; ++i) {
        const std::string arg = argv[i];
        auto nextDouble = [&](const char* flag) -> double {
            if (i + 1 >= argc)
                throw std::runtime_error(std::string(flag) + " requires a value.");
            return std::stod(argv[++i]);
        };
        auto nextInt = [&](const char* flag) -> int {
            if (i + 1 >= argc)
                throw std::runtime_error(std::string(flag) + " requires a value.");
            return std::stoi(argv[++i]);
        };
        if      (arg == "--dos")  doDos  = true;
        else if (arg == "--band") doBand = true;
        else if (arg == "--ha")   doHa   = true;
        else if (arg == "--tmin")  tMin  = nextDouble("--tmin");
        else if (arg == "--tmax")  tMax  = nextDouble("--tmax");
        else if (arg == "--dt")    dt    = nextDouble("--dt");
        else if (arg == "--natom") natom = nextInt("--natom");
        else if (arg == "--labels") {
            if (i + 1 >= argc) throw std::runtime_error("--labels requires a list.");
            manualLabels = parse_labels(argv[++i]);
        } else if (arg == "--qlabels") {
            if (i + 1 >= argc) throw std::runtime_error("--qlabels requires a file path.");
            qLabelsFile = argv[++i];
        } else if (!outPrefixSet && arg.rfind("--", 0) != 0) {
            outPrefix = arg;
            outPrefixSet = true;
        } else {
            throw std::runtime_error("Unknown argument for 'phonon -post': " + arg);
        }
    }

    // Default: run all three if none explicitly selected
    if (!doDos && !doBand && !doHa) { doDos = doBand = doHa = true; }

    if (natom < 1) throw std::runtime_error("--natom must be >= 1.");
    if (tMax <= tMin || dt <= 0.0)
        throw std::runtime_error("Invalid temperature range.");

    const std::string dosFile  = prefix + ".phonon.dos";
    const std::string bandFile = prefix + ".phonon_band.freq";

    if (doDos || doHa) {
        const auto dosData = parse_matdyn_dos(dosFile);
        if (doDos)  write_phonon_dos_plot(dosData, outPrefix);
        if (doHa) {
            const auto ha = compute_harmonic_approximation(
                dosData, make_temp_grid(tMin, tMax, dt), natom);
            write_ha_report(ha, outPrefix);
            write_ha_plot(ha, outPrefix);
        }
    }

    if (doBand) {
        auto bandData = parse_matdyn_freq(bandFile, qLabelsFile);
        if (!manualLabels.empty()) {
            if (manualLabels.size() == bandData.labelMarks.size()) {
                for (size_t k = 0; k < manualLabels.size(); ++k)
                    bandData.labelMarks[k].second = manualLabels[k];
            } else {
                // Size mismatch: assign first/last and evenly-spaced interior labels
                std::cerr << "  WARNING: " << manualLabels.size() << " labels given but "
                          << bandData.labelMarks.size() << " corners detected; "
                          << "assigning labels to first/last and interior corners.\n";
                // Rebuild labelMarks to match manualLabels count
                const int nLabels = static_cast<int>(manualLabels.size());
                const int nMark   = static_cast<int>(bandData.labelMarks.size());
                std::vector<std::pair<double, std::string>> newMarks;
                newMarks.push_back({bandData.distances.front(), manualLabels[0]});
                // Distribute intermediate labels across intermediate marks
                for (int li = 1; li < nLabels - 1; ++li) {
                    int mi = (nMark > 2)
                        ? std::clamp(1 + (li - 1) * (nMark - 2) / std::max(nLabels - 2, 1),
                                     1, nMark - 2)
                        : 0;
                    if (mi > 0 && mi < nMark)
                        newMarks.push_back({bandData.labelMarks[mi].first, manualLabels[li]});
                    else
                        newMarks.push_back({bandData.distances.back() * li / (nLabels - 1),
                                            manualLabels[li]});
                }
                newMarks.push_back({bandData.distances.back(), manualLabels[nLabels - 1]});
                bandData.labelMarks = std::move(newMarks);
            }
        }
        write_phonon_band_plot(bandData, outPrefix);
    }

    return 0;
}

}  // namespace qe
