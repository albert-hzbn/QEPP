#include "qe/cli/commands.hpp"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "qe/band.hpp"
#include "qe/cif.hpp"
#include "qe/dos.hpp"
#include "qe/help.hpp"
#include "qe/qe_input.hpp"
#include "qe/types.hpp"
#include "qe/utils.hpp"

namespace qe {

int handle_band_post_mode(int argc, char** argv, int s) {
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

int handle_band_fat_mode(int argc, char** argv, int s) {
    // Usage: band -fat <bands.gnu> <atomic_proj.xml> [outprefix] [fermi] [filter ...]
    // filter: element=Si,Fe  |  atom=1,2,3  |  orbital=s,p,d,f  (combinable)
    if (argc < 4 + s) {
        print_help_command(argv[0], "band", "-fat");
        return 1;
    }

    const std::string bandPath = argv[2 + s];
    const std::string xmlPath = argv[3 + s];
    std::string outPrefix = stem_from_path(bandPath);
    double fermiEv = 0.0;
    std::vector<std::string> elements;
    std::vector<int> atomNums;
    std::vector<int> orbitalLs;

    auto l_from_name = [](const std::string& name) -> int {
        const std::string n = qe::to_lower(name);
        if (n == "s") return 0;
        if (n == "p") return 1;
        if (n == "d") return 2;
        if (n == "f") return 3;
        return -1;
    };

    // Optional positional args: [outprefix] [fermi] [filter ...]
    for (int i = 4 + s; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg.rfind("element=", 0) == 0) {
            std::istringstream ss(arg.substr(8));
            std::string tok;
            while (std::getline(ss, tok, ','))
                if (!tok.empty()) elements.push_back(trim(tok));
        } else if (arg.rfind("atom=", 0) == 0) {
            std::istringstream ss(arg.substr(5));
            std::string tok;
            while (std::getline(ss, tok, ','))
                if (!tok.empty()) atomNums.push_back(std::stoi(trim(tok)));
        } else if (arg.rfind("orbital=", 0) == 0) {
            std::istringstream ss(arg.substr(8));
            std::string tok;
            while (std::getline(ss, tok, ',')) {
                const int l = l_from_name(trim(tok));
                if (l >= 0) orbitalLs.push_back(l);
            }
        } else {
            double v = 0.0;
            if (try_parse_double(arg, v)) {
                fermiEv = v;
            } else if (i == 4 + s) {
                // Could be fermi from QE output file or outprefix
                // Try as QE output for Fermi energy extraction first
                try {
                    fermiEv = extract_fermi_from_qe_output(arg);
                    std::cout << "Extracted Fermi energy from QE output: "
                              << std::fixed << std::setprecision(6) << fermiEv << " eV\n";
                } catch (...) {
                    outPrefix = arg;
                }
            } else {
                outPrefix = arg;
            }
        }
    }

    auto bandData = parse_band_table(bandPath);
    const auto proj = parse_atomic_proj(xmlPath);
    const auto groups = build_fatband_groups(proj, elements, atomNums, orbitalLs);
    write_fatband_plots(bandData, groups, fermiEv, outPrefix);
    return 0;
}

int handle_band_pre_mode(int argc, char** argv, int s) {
    if (argc < 4 + s || argc > 8 + s) {
        print_help_command(argv[0], "band", "-pre");
        return 1;
    }

    const std::string cifPath = argv[2 + s];
    const std::string scfInputPath = argv[3 + s];
    const std::string defaultPrefix = stem_from_path(cifPath);
    const std::string bandsPwPath = (argc >= 5 + s) ? argv[4 + s] : (defaultPrefix + ".bands.in");
    const std::string bandsPpPath = (argc >= 6 + s) ? argv[5 + s] : (defaultPrefix + ".bands_pp.in");
    const int pointsPerSegment = (argc >= 7 + s) ? std::stoi(argv[6 + s]) : 20;
    const int nbnd = (argc >= 8 + s) ? std::stoi(argv[7 + s]) : 0;

    const CifStructure structure = parse_cif(cifPath);
    const SymmetryKPath kpath = suggest_kpath_from_cif(structure);

    write_bands_input_from_scf_template(scfInputPath, bandsPwPath, kpath,
                                        pointsPerSegment, nbnd);
    const std::string filbandName = defaultPrefix + ".bands.dat";
    write_bands_pp_input_from_scf_template(scfInputPath, bandsPpPath, filbandName);

    // Also write a projwfc.x input alongside bands_pp for fat-band post-processing
    const std::string projwfcPath = stem_from_path(bandsPpPath) + ".projwfc.in";
    write_bands_projwfc_input(scfInputPath, projwfcPath);

    std::cout << "Generated QE bands pw.x input: " << bandsPwPath << "\n";
    std::cout << "Generated QE bands.x input: " << bandsPpPath << "\n";
    std::cout << "Generated projwfc.x input: " << projwfcPath << "\n";
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

}  // namespace qe
