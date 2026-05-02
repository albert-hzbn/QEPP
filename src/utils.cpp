#include "qe/utils.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <map>
#include <stdexcept>

namespace qe {

std::string trim(const std::string& s) {
    const auto start = s.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) {
        return "";
    }
    const auto end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}

std::string strip_quotes(const std::string& s) {
    if (s.size() >= 2 &&
        ((s.front() == '\'' && s.back() == '\'') ||
         (s.front() == '"' && s.back() == '"'))) {
        return s.substr(1, s.size() - 2);
    }
    return s;
}

std::vector<std::string> split_cif_row(const std::string& line) {
    std::vector<std::string> tokens;
    std::string current;
    bool inQuote = false;
    char quoteChar = '\0';

    for (char c : line) {
        if (!inQuote && (c == '\'' || c == '"')) {
            inQuote = true;
            quoteChar = c;
            current.push_back(c);
        } else if (inQuote && c == quoteChar) {
            inQuote = false;
            current.push_back(c);
        } else if (!inQuote && std::isspace(static_cast<unsigned char>(c))) {
            if (!current.empty()) {
                tokens.push_back(strip_quotes(current));
                current.clear();
            }
        } else {
            current.push_back(c);
        }
    }
    if (!current.empty()) {
        tokens.push_back(strip_quotes(current));
    }

    return tokens;
}

double parse_double(const std::string& token) {
    const auto parenPos = token.find('(');
    const std::string cleaned =
        (parenPos == std::string::npos) ? token : token.substr(0, parenPos);
    return std::stod(cleaned);
}

std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

std::string normalize_symbol(const std::string& raw) {
    if (raw.empty()) {
        return raw;
    }

    std::string letters;
    for (char c : raw) {
        if (std::isalpha(static_cast<unsigned char>(c))) {
            letters.push_back(c);
        } else {
            break;
        }
    }

    if (letters.empty()) {
        return raw;
    }

    std::string out;
    out.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(letters[0]))));
    if (letters.size() > 1) {
        out.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(letters[1]))));
    }
    return out;
}

std::string stem_from_path(const std::string& path) {
    const size_t slash = path.find_last_of("/");
    const std::string file = (slash == std::string::npos) ? path : path.substr(slash + 1);
    const size_t dot = file.find_last_of('.');
    return (dot == std::string::npos) ? file : file.substr(0, dot);
}

std::vector<std::string> load_lines(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open()) {
        throw std::runtime_error("Could not open file: " + path);
    }

    std::vector<std::string> out;
    std::string line;
    while (std::getline(in, line)) {
        out.push_back(line);
    }
    return out;
}

bool try_parse_double(const std::string& s, double& value) {
    try {
        size_t consumed = 0;
        value = std::stod(s, &consumed);
        return consumed == s.size();
    } catch (...) {
        return false;
    }
}

std::string extract_quoted_assignment(const std::vector<std::string>& lines,
                                      const std::string& keyLower) {
    for (const auto& raw : lines) {
        const std::string lower = to_lower(raw);
        const auto keyPos = lower.find(keyLower);
        if (keyPos == std::string::npos) {
            continue;
        }

        const auto eqPos = raw.find('=', keyPos);
        if (eqPos == std::string::npos) {
            continue;
        }

        const auto q1 = raw.find('"', eqPos);
        if (q1 != std::string::npos) {
            const auto q2 = raw.find('"', q1 + 1);
            if (q2 != std::string::npos) {
                return raw.substr(q1 + 1, q2 - q1 - 1);
            }
        }

        const auto s1 = raw.find('\'', eqPos);
        if (s1 != std::string::npos) {
            const auto s2 = raw.find('\'', s1 + 1);
            if (s2 != std::string::npos) {
                return raw.substr(s1 + 1, s2 - s1 - 1);
            }
        }
    }

    return "";
}

double atomic_mass(const std::string& symbol) {
    // Atomic masses in amu (IUPAC 2021 standard values)
    static const std::map<std::string, double> masses = {
        {"H",1.008},   {"He",4.0026},  {"Li",6.94},    {"Be",9.0122},
        {"B",10.81},   {"C",12.011},   {"N",14.007},   {"O",15.999},
        {"F",18.998},  {"Ne",20.180},  {"Na",22.990},  {"Mg",24.305},
        {"Al",26.982}, {"Si",28.085},  {"P",30.974},   {"S",32.06},
        {"Cl",35.45},  {"Ar",39.948}, {"K",39.098},   {"Ca",40.078},
        {"Sc",44.956}, {"Ti",47.867}, {"V",50.942},   {"Cr",51.996},
        {"Mn",54.938}, {"Fe",55.845}, {"Co",58.933},  {"Ni",58.693},
        {"Cu",63.546}, {"Zn",65.38},  {"Ga",69.723},  {"Ge",72.630},
        {"As",74.922}, {"Se",78.971}, {"Br",79.904},  {"Kr",83.798},
        {"Rb",85.468}, {"Sr",87.62},  {"Y",88.906},   {"Zr",91.224},
        {"Nb",92.906}, {"Mo",95.95},  {"Tc",97.0},    {"Ru",101.07},
        {"Rh",102.91}, {"Pd",106.42}, {"Ag",107.87},  {"Cd",112.41},
        {"In",114.82}, {"Sn",118.71}, {"Sb",121.76},  {"Te",127.60},
        {"I",126.90},  {"Xe",131.29}, {"Cs",132.91},  {"Ba",137.33},
        {"La",138.91}, {"Ce",140.12}, {"Pr",140.91},  {"Nd",144.24},
        {"Pm",145.0},  {"Sm",150.36}, {"Eu",151.96},  {"Gd",157.25},
        {"Tb",158.93}, {"Dy",162.50}, {"Ho",164.93},  {"Er",167.26},
        {"Tm",168.93}, {"Yb",173.05}, {"Lu",174.97},  {"Hf",178.49},
        {"Ta",180.95}, {"W",183.84},  {"Re",186.21},  {"Os",190.23},
        {"Ir",192.22}, {"Pt",195.08}, {"Au",196.97},  {"Hg",200.59},
        {"Tl",204.38}, {"Pb",207.2},  {"Bi",208.98},  {"Po",209.0},
        {"At",210.0},  {"Rn",222.0},  {"Fr",223.0},   {"Ra",226.0},
        {"Ac",227.0},  {"Th",232.04}, {"Pa",231.04},  {"U",238.03},
        {"Np",237.0},  {"Pu",244.0},  {"Am",243.0},   {"Cm",247.0},
        {"Bk",247.0},  {"Cf",251.0},  {"Es",252.0},   {"Fm",257.0},
        {"Md",258.0},  {"No",259.0},  {"Lr",262.0}
    };
    const auto it = masses.find(symbol);
    return (it != masses.end()) ? it->second : 28.085;  // default to Si mass
}

}  // namespace qe
