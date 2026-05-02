#include "qe/cli/commands.hpp"

#include <string>

#include "qe/elastic.hpp"
#include "qe/help.hpp"

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
