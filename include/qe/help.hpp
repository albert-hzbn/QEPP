#pragma once

#include <string>

void print_help(const char* prog);
void print_help_command(const char* prog, const std::string& cmd,
                        const std::string& sub = "");
