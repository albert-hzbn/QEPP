#pragma once

#include <string>
#include <vector>

namespace qe {

std::string trim(const std::string& s);
std::string strip_quotes(const std::string& s);
std::vector<std::string> split_cif_row(const std::string& line);
double parse_double(const std::string& token);
std::string to_lower(std::string s);
std::string normalize_symbol(const std::string& raw);
std::string stem_from_path(const std::string& path);
std::vector<std::string> split_csv(const std::string& text);
std::vector<std::string> load_lines(const std::string& path);
bool try_parse_double(const std::string& s, double& value);
std::string extract_quoted_assignment(const std::vector<std::string>& lines,
                                      const std::string& keyLower);
bool is_directory(const std::string& path);
std::string join_paths(const std::string& base, const std::string& leaf);

// Return atomic mass in amu for an element symbol; returns 28.085 (Si) if unknown.
double atomic_mass(const std::string& symbol);

}  // namespace qe
