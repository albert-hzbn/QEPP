#include "qe/cli/commands.hpp"

#include <string>

#include "qe/elastic.hpp"
#include "qe/help.hpp"

namespace {

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
    if (argc < 4 + s) {
        print_help_command(argv[0], "elastic", "-pre");
        return 1;
    }
    const std::string tmpl = argv[2 + s];
    const std::string outDir = argv[3 + s];
    const int nDeltas = (argc >= 5 + s) ? std::stoi(argv[4 + s]) : 7;
    const double maxDelta = (argc >= 6 + s) ? std::stod(argv[5 + s]) : 0.04;
    generate_elastic_inputs(tmpl, outDir, nDeltas, maxDelta);
    return 0;
}

int handle_elastic_run_mode(int argc, char** argv, int s) {
    if (argc < 4 + s) {
        print_help_command(argv[0], "elastic", "-run");
        return 1;
    }
    const std::string tmpl = argv[2 + s];
    const std::string outDir = argv[3 + s];
    const QeParallelOptions parallel = parse_parallel_options(argc, argv, 4 + s);
    run_elastic_dataset(tmpl, outDir, parallel);
    return 0;
}

int handle_elastic_post_mode(int argc, char** argv, int s) {
    if (argc < 4 + s) {
        print_help_command(argv[0], "elastic", "-post");
        return 1;
    }
    const std::string tmpl = argv[2 + s];
    const std::string outDir = argv[3 + s];
    const ElasticResults res = compute_elastic_properties(tmpl, outDir);
    print_elastic_report(res, res.crystalFamily, outDir + "/elastic_results.txt");
    return 0;
}

}  // namespace qe
