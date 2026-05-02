#include "qe/band.hpp"

#include <matplot/matplot.h>

#include <cctype>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

#include "qe/utils.hpp"

namespace qe {

static bool extract_header_int_after_key(const std::string& line,
                                         const std::string& key,
                                         int& value) {
    const std::string lower = to_lower(line);
    const auto pos = lower.find(key);
    if (pos == std::string::npos) {
        return false;
    }

    size_t i = pos + key.size();
    while (i < line.size() && (line[i] == ' ' || line[i] == '\t' || line[i] == '=')) {
        ++i;
    }

    size_t j = i;
    while (j < line.size() && std::isdigit(static_cast<unsigned char>(line[j]))) {
        ++j;
    }

    if (j == i) {
        return false;
    }

    value = std::stoi(line.substr(i, j - i));
    return true;
}

BandData parse_band_table(const std::string& bandPath) {
    std::ifstream in(bandPath);
    if (!in.is_open()) {
        throw std::runtime_error("Could not open band data file: " + bandPath);
    }

    std::vector<std::string> lines;
    std::string line;
    while (std::getline(in, line)) {
        lines.push_back(line);
    }

    int nbnd = -1;
    int nks = -1;
    size_t headerIdx = 0;
    for (size_t i = 0; i < lines.size(); ++i) {
        int a = -1;
        int b = -1;
        const bool hasNbnd = extract_header_int_after_key(lines[i], "nbnd", a);
        const bool hasNks = extract_header_int_after_key(lines[i], "nks", b);
        if (hasNbnd && hasNks) {
            nbnd = a;
            nks = b;
            headerIdx = i;
            break;
        }
    }

    if (nbnd > 0 && nks > 0) {
        std::vector<Eigen::Vector3d> kvec;
        kvec.reserve(static_cast<size_t>(nks));
        std::vector<std::vector<double>> energies(
            static_cast<size_t>(nbnd), std::vector<double>(static_cast<size_t>(nks), 0.0));

        size_t i = headerIdx + 1;
        for (int ik = 0; ik < nks; ++ik) {
            std::vector<double> coords;
            while (i < lines.size()) {
                const std::string t = trim(lines[i++]);
                if (t.empty()) {
                    continue;
                }
                std::istringstream iss(t);
                double v = 0.0;
                coords.clear();
                while (iss >> v) {
                    coords.push_back(v);
                }
                if (coords.size() >= 3) {
                    break;
                }
            }

            if (coords.size() < 3) {
                throw std::runtime_error(
                    "Unexpected end of band file while reading k-point coordinates.");
            }

            kvec.emplace_back(coords[0], coords[1], coords[2]);

            std::vector<double> e;
            e.reserve(static_cast<size_t>(nbnd));
            while (i < lines.size() && static_cast<int>(e.size()) < nbnd) {
                const std::string t = trim(lines[i++]);
                if (t.empty()) {
                    continue;
                }
                std::istringstream iss(t);
                double v = 0.0;
                while (iss >> v) {
                    e.push_back(v);
                    if (static_cast<int>(e.size()) == nbnd) {
                        break;
                    }
                }
            }

            if (static_cast<int>(e.size()) < nbnd) {
                throw std::runtime_error(
                    "Unexpected end of band file while reading energies.");
            }

            for (int b = 0; b < nbnd; ++b) {
                energies[static_cast<size_t>(b)][static_cast<size_t>(ik)] =
                    e[static_cast<size_t>(b)];
            }
        }

        std::vector<double> kdist(static_cast<size_t>(nks), 0.0);
        for (int ik = 1; ik < nks; ++ik) {
            const double dk =
                (kvec[static_cast<size_t>(ik)] - kvec[static_cast<size_t>(ik - 1)]).norm();
            kdist[static_cast<size_t>(ik)] = kdist[static_cast<size_t>(ik - 1)] + dk;
        }

        // Detect high-symmetry segment boundaries by direction changes.
        std::vector<std::pair<double, std::string>> marks;
        marks.emplace_back(kdist[0], "");
        for (int ik = 1; ik < nks - 1; ++ik) {
            const Eigen::Vector3d dkPrev =
                kvec[static_cast<size_t>(ik)] - kvec[static_cast<size_t>(ik - 1)];
            const Eigen::Vector3d dkNext =
                kvec[static_cast<size_t>(ik + 1)] - kvec[static_cast<size_t>(ik)];
            const double normPrev = dkPrev.norm();
            const double normNext = dkNext.norm();
            if (normPrev > 1e-10 && normNext > 1e-10) {
                const double cosAngle = dkPrev.dot(dkNext) / (normPrev * normNext);
                if (cosAngle < 0.999) {
                    marks.emplace_back(kdist[static_cast<size_t>(ik)], "");
                }
            }
        }
        marks.emplace_back(kdist[static_cast<size_t>(nks - 1)], "");

        BandData out;
        out.kByBand.resize(static_cast<size_t>(nbnd), kdist);
        out.eByBand = std::move(energies);
        out.kLabelMarks = std::move(marks);
        return out;
    }

    BandData out;
    std::vector<double> kCurrent;
    std::vector<double> eCurrent;

    auto flush_current = [&]() {
        if (!kCurrent.empty()) {
            out.kByBand.push_back(kCurrent);
            out.eByBand.push_back(eCurrent);
            kCurrent.clear();
            eCurrent.clear();
        }
    };

    for (const auto& raw : lines) {
        const std::string t = trim(raw);
        if (t.empty()) {
            flush_current();
            continue;
        }

        if (t[0] == '#' || t[0] == '!' || t[0] == '@') {
            continue;
        }

        std::istringstream iss(t);
        double k = 0.0;
        double e = 0.0;
        if (iss >> k >> e) {
            kCurrent.push_back(k);
            eCurrent.push_back(e);
        }
    }

    flush_current();

    if (out.kByBand.empty()) {
        throw std::runtime_error(
            "No band data found. Expected a 2-column file (k, E) with blank lines "
            "between bands.");
    }

    return out;
}

void write_band_plot_bundle(const std::string& bandInputPath,
                            const BandData& bandData,
                            double fermiEv,
                            const std::string& outPrefix) {
    const std::string dataPath = outPrefix + ".band.dat";
    const std::string pngPath = outPrefix + ".band.png";

    std::ofstream dataOut(dataPath);
    if (!dataOut.is_open()) {
        throw std::runtime_error("Could not open band output data file: " + dataPath);
    }

    dataOut << std::fixed << std::setprecision(10);
    for (size_t b = 0; b < bandData.kByBand.size(); ++b) {
        const auto& k = bandData.kByBand[b];
        const auto& e = bandData.eByBand[b];
        for (size_t i = 0; i < k.size(); ++i) {
            dataOut << k[i] << " " << (e[i] - fermiEv) << "\n";
        }
        dataOut << "\n";
    }
    dataOut.close();

    using namespace matplot;
    auto fig = figure(true);
    fig->size(1400, 900);
    hold(on);

    double kMin = std::numeric_limits<double>::infinity();
    double kMax = -std::numeric_limits<double>::infinity();
    double eMin = std::numeric_limits<double>::infinity();
    double eMax = -std::numeric_limits<double>::infinity();

    // First pass: compute k and energy ranges from actual data.
    for (size_t b = 0; b < bandData.kByBand.size(); ++b) {
        const auto& k = bandData.kByBand[b];
        for (size_t i = 0; i < k.size(); ++i) {
            kMin = std::min(kMin, k[i]);
            kMax = std::max(kMax, k[i]);
            const double e = bandData.eByBand[b][i] - fermiEv;
            eMin = std::min(eMin, e);
            eMax = std::max(eMax, e);
        }
    }

    // Add padding and clip to a sensible window (±15 eV around Fermi).
    const double yPad = (eMax - eMin) * 0.05 + 0.5;
    const double yLo = std::max(eMin - yPad, -15.0);
    const double yHi = std::min(eMax + yPad,  15.0);

    // Draw vertical lines first so bands render on top.
    if (!bandData.kLabelMarks.empty()) {
        for (const auto& [kPos, label] : bandData.kLabelMarks) {
            if (kPos > kMin + 1e-10 && kPos < kMax - 1e-10) {
                auto vLine = plot(std::vector<double>{kPos, kPos},
                                  std::vector<double>{yLo, yHi}, "--");
                // color_array is RGBA — must include alpha=1.
                vLine->color({0.55F, 0.55F, 0.55F, 1.0F});
                vLine->line_width(1.0);
            }
        }
    }

    // Fermi level line.
    auto efLine = plot(std::vector<double>{kMin, kMax}, std::vector<double>{0.0, 0.0}, "--");
    efLine->color({0.2F, 0.2F, 0.2F, 1.0F});
    efLine->line_width(1.0);

    // Band curves.
    for (size_t b = 0; b < bandData.kByBand.size(); ++b) {
        const auto& k = bandData.kByBand[b];
        std::vector<double> eShifted;
        eShifted.reserve(bandData.eByBand[b].size());
        for (size_t i = 0; i < k.size(); ++i) {
            eShifted.push_back(bandData.eByBand[b][i] - fermiEv);
        }
        auto p = plot(k, eShifted);
        p->line_width(1.8);
        // color_array is RGBA — alpha must be 1.0 or the line is transparent.
        p->color({0.10F, 0.24F, 0.56F, 1.0F});
    }

    // X-axis: high-symmetry labels.
    if (!bandData.kLabelMarks.empty()) {
        std::vector<double> tickPositions;
        std::vector<std::string> tickLabels;
        for (const auto& [kPos, label] : bandData.kLabelMarks) {
            tickPositions.push_back(kPos);
            std::string displayLabel = label;
            if (displayLabel == "G" || displayLabel == "GM" || displayLabel == "Gamma") {
                displayLabel = "{/Symbol G}";
            }
            tickLabels.push_back(displayLabel);
        }
        xticks(tickPositions);
        xticklabels(tickLabels);
    }

    xlim({kMin, kMax});
    ylim({yLo, yHi});

    title("Quantum ESPRESSO Band Structure");
    xlabel("");
    ylabel("E - E_{F} (eV)");
    grid(on);

    std::cout << "Prepared band data: " << dataPath << "\n";
    save(pngPath);
    std::cout << "Generated band plot image (Matplot++): " << pngPath << "\n";
    std::cout << "Source band file: " << bandInputPath << "\n";
}

}  // namespace qe
