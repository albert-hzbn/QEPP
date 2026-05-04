#include "qe/cli/commands.hpp"

#include <cctype>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

#include "qe/elastic.hpp"
#include "qe/help.hpp"
#include "qe/utils.hpp"

namespace {

namespace fs = std::filesystem;

bool is_parallel_flag(const std::string& arg) {
    return arg == "--np" || arg == "--ni" || arg == "--nk" ||
           arg == "--nb" || arg == "--nt" || arg == "--nd";
}

std::string detect_template_from_outdir(const std::string& outDir) {
    const fs::path out(outDir);
    const fs::path parent = out.parent_path().empty() ? fs::path(".") : out.parent_path();

    std::vector<fs::path> parentInputs;
    if (fs::is_directory(parent)) {
        for (const auto& e : fs::directory_iterator(parent)) {
            if (e.is_regular_file() && e.path().extension() == ".in")
                parentInputs.push_back(e.path());
        }
    }
    if (parentInputs.size() == 1)
        return parentInputs.front().string();

    std::vector<fs::path> cwdInputs;
    for (const auto& e : fs::directory_iterator(fs::path("."))) {
        if (e.is_regular_file() && e.path().extension() == ".in")
            cwdInputs.push_back(e.path());
    }
    if (cwdInputs.size() == 1)
        return cwdInputs.front().string();

    throw std::runtime_error(
        "Could not auto-detect scf_template.in for elastic folder '" + outDir +
        "'. Provide explicit arguments: elastic <mode> <scf_template.in> <outdir>");
}

std::pair<std::string, std::string> resolve_elastic_paths(const std::vector<std::string>& positional) {
    if (positional.size() >= 2)
        return {positional[0], positional[1]};

    fs::path base = positional.empty() ? fs::path(".") : fs::path(positional[0]);
    if (!fs::is_directory(base))
        throw std::runtime_error("Expected an elastic directory or volume directory: " + base.string());

    fs::path outDir = base;
    if (fs::is_directory(base / "elastic"))
        outDir = base / "elastic";

    const std::string tmpl = detect_template_from_outdir(outDir.string());
    return {tmpl, outDir.string()};
}

qe::QeParallelOptions parse_parallel_options(int argc, char** argv, int start) {
    qe::QeParallelOptions opts;
    for (int i = start; i < argc; ++i) {
        const std::string arg = argv[i];
        auto nextInt = [&](const std::string& flag) -> int {
            if (i + 1 >= argc)
                throw std::runtime_error(flag + " requires an integer value.");
            return std::stoi(argv[++i]);
        };
        if      (arg == "--np") opts.np = nextInt(arg);
        else if (arg == "--ni") opts.ni = nextInt(arg);
        else if (arg == "--nk") opts.nk = nextInt(arg);
        else if (arg == "--nb") opts.nb = nextInt(arg);
        else if (arg == "--nt") opts.nt = nextInt(arg);
        else if (arg == "--nd") opts.nd = nextInt(arg);
        else throw std::runtime_error("Unknown argument for elastic run: " + arg);
    }
    return opts;
}

}  // namespace

namespace qe {

int handle_elastic_pre_mode(int argc, char** argv, int s) {
    if (argc < 3 + s) {
        print_help_command(argv[0], "elastic", "-pre");
        return 1;
    }

    const std::string tmpl = argv[2 + s];
    std::string outDir;
    int nDeltas = 7;
    double maxDelta = 0.04;

    int i = 3 + s;
    if (i < argc && std::string(argv[i]).rfind("--", 0) != 0) {
        outDir = argv[i++];
    }
    if (i < argc && std::string(argv[i]).rfind("--", 0) != 0) {
        nDeltas = std::stoi(argv[i++]);
    }
    if (i < argc && std::string(argv[i]).rfind("--", 0) != 0) {
        maxDelta = std::stod(argv[i++]);
    }

    for (; i < argc; ++i) {
        const std::string arg = argv[i];
        auto nextInt = [&](const std::string& flag) -> int {
            if (i + 1 >= argc)
                throw std::runtime_error(flag + " requires an integer value.");
            return std::stoi(argv[++i]);
        };
        auto nextDouble = [&](const std::string& flag) -> double {
            if (i + 1 >= argc)
                throw std::runtime_error(flag + " requires a numeric value.");
            return std::stod(argv[++i]);
        };
        auto nextString = [&](const std::string& flag) -> std::string {
            if (i + 1 >= argc)
                throw std::runtime_error(flag + " requires a value.");
            return argv[++i];
        };

        if      (arg == "--outdir")   outDir   = nextString(arg);
        else if (arg == "--ndeltas")  nDeltas  = nextInt(arg);
        else if (arg == "--maxdelta") maxDelta = nextDouble(arg);
        else throw std::runtime_error("Unknown argument for 'elastic -pre': " + arg);
    }

    if (outDir.empty())
        outDir = stem_from_path(tmpl) + "_elastic";

    generate_elastic_inputs(tmpl, outDir, nDeltas, maxDelta);
    return 0;
}

int handle_elastic_run_mode(int argc, char** argv, int s) {
    if (argc < 3 + s) {
        print_help_command(argv[0], "elastic", "-run");
        return 1;
    }

    std::vector<std::string> positional;
    int i = 2 + s;
    while (i < argc && std::string(argv[i]).rfind("--", 0) != 0) {
        positional.emplace_back(argv[i]);
        ++i;
    }
    for (const auto& p : positional) {
        if (is_parallel_flag(p))
            throw std::runtime_error("Invalid positional argument: " + p);
    }

    const auto paths = resolve_elastic_paths(positional);
    const std::string& tmpl = paths.first;
    const std::string& outDir = paths.second;
    const QeParallelOptions parallel = parse_parallel_options(argc, argv, i);
    run_elastic_dataset(tmpl, outDir, parallel);
    return 0;
}

int handle_elastic_post_mode(int argc, char** argv, int s) {
    if (argc < 3 + s) {
        print_help_command(argv[0], "elastic", "-post");
        return 1;
    }

    std::vector<std::string> positional;
    for (int i = 2 + s; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg.rfind("--", 0) == 0)
            throw std::runtime_error("Unknown option for 'elastic -post': " + arg);
        positional.push_back(arg);
    }

    const auto paths = resolve_elastic_paths(positional);
    const std::string& tmpl = paths.first;
    const std::string& outDir = paths.second;
    const ElasticResults res = compute_elastic_properties(tmpl, outDir);
    print_elastic_report(res, res.crystalFamily, outDir + "/elastic_results.txt");
    return 0;
}

}  // namespace qe
