#include <iostream>
#include <string>
#include <windows.h>

using namespace std;
#define PRINTED_FIELD_SIZE 4

#ifndef PEA_0_MAIN_PROGRAM_H
#define PEA_0_MAIN_PROGRAM_H

string get_current_dir() {
    char result[MAX_PATH];
    return string(result, GetModuleFileName(nullptr, result, MAX_PATH));
}

string current_dir = get_current_dir();
string tsp_path = R"(\..\..\..\Data\TSP\)";
string small_path = R"(\..\..\..\Data\SMALL\)";
string atsp_path = R"(\..\..\..\Data\ATSP\)";

#endif //PEA_0_MAIN_PROGRAM_H
