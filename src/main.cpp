#include <iostream>
#include <stdexcept>

#include "qe/cli/dispatch.hpp"

int main(int argc, char** argv) {
    try {
        return qe::dispatch_cli(argc, argv);
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << "\n";
        return 1;
    }
}
