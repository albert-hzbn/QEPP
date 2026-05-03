#include "qe/cli/commands.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

#include "qe/help.hpp"
#include "qe/struct.hpp"
#include "qe/utils.hpp"

namespace qe {

int handle_struct_mode(int argc, char** argv, int s) {
    // struct -post <structure_file> [output_prefix] [--sro] [--nshells N]
    //               [--tol T] [--source input|output|auto]
    if (argc < 3 + s) {
        print_help_command(argv[0], "struct", "-post");
        return 1;
    }

    const std::string structurePath = argv[2 + s];
    std::string outPrefix = stem_from_path(structurePath);

    bool doSro = false;
    int nShells = 2;
    double shellTol = 1e-3;
    std::string source = "auto";

    bool outPrefixSet = false;
    for (int i = 3 + s; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--sro") {
            doSro = true;
            continue;
        }
        if (arg == "--nshells") {
            if (i + 1 >= argc)
                throw std::runtime_error("--nshells requires an integer value.");
            nShells = std::stoi(argv[++i]);
            doSro = true;
            continue;
        }
        if (arg == "--tol") {
            if (i + 1 >= argc)
                throw std::runtime_error("--tol requires a numeric value.");
            shellTol = std::stod(argv[++i]);
            doSro = true;
            continue;
        }
        if (arg == "--source") {
            if (i + 1 >= argc)
                throw std::runtime_error("--source requires one of: input, output, auto.");
            source = to_lower(argv[++i]);
            if (source != "input" && source != "output" && source != "auto") {
                throw std::runtime_error("Invalid --source value: '" + source +
                                         "'. Use input, output, or auto.");
            }
            continue;
        }

        if (!outPrefixSet && arg.rfind("--", 0) != 0) {
            outPrefix = arg;
            outPrefixSet = true;
            continue;
        }

        throw std::runtime_error("Unknown argument for struct -post: " + arg);
    }

    auto valid = [](const StructInfo& info) {
        return !info.atoms.empty() && std::abs(info.cellAngst.determinant()) > 1e-12;
    };

    StructInfo info;
    if (source == "input") {
        info = parse_struct_from_qe_input(structurePath);
    } else if (source == "output") {
        info = parse_struct_from_qe_output(structurePath);
    } else {
        const std::string lowPath = to_lower(structurePath);
        const bool looksOutput = lowPath.find(".out") != std::string::npos;
        if (looksOutput) {
            info = parse_struct_from_qe_output(structurePath);
            if (!valid(info)) info = parse_struct_from_qe_input(structurePath);
        } else {
            info = parse_struct_from_qe_input(structurePath);
            if (!valid(info)) info = parse_struct_from_qe_output(structurePath);
        }
    }

    if (!valid(info))
        throw std::runtime_error(
            "Could not parse structure from: " + structurePath +
            "\n  Provide a QE SCF input with ATOMIC_POSITIONS/CELL_PARAMETERS, or "
            "a QE output containing complete final CELL_PARAMETERS/ATOMIC_POSITIONS blocks.");

    write_struct_report(info, outPrefix);

    if (doSro) {
        const SroReport report = estimate_warren_cowley_sro(info, nShells, shellTol);
        write_sro_report(report, outPrefix);
        const SroReportRC rcReport = estimate_sro_rao_curtin(info, nShells, shellTol);
        write_sro_rc_report(rcReport, outPrefix);
    }
    return 0;
}

}  // namespace qe
