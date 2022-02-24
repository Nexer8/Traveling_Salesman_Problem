#include <iostream>
#include <string>
#include <vector>

using namespace std;

#ifndef PEA_0_MAIN_PROGRAM_H
#define PEA_0_MAIN_PROGRAM_H

#define PRINTED_FIELD_SIZE 4

enum BeginningSolution {
    NATURAL,
    RANDOM,
    GREEDY
};

enum NeighborhoodType {
    SWAP,
    INVERT,
    INSERT
};

enum ProblemType {
    TSP,
    SMALL,
    ATSP
};

typedef int (*move_linear_update)(int &, vector<int> &, int &, int &, vector<vector<int>> &);

typedef void (*move_foo)(vector<int> &, int, int);

constexpr string_view ATSP_PATH("../data/ATSP/");
constexpr string_view SMALL_PATH = "../data/SMALL/";
constexpr string_view TSP_PATH = "../data/TSP/";

#endif //PEA_0_MAIN_PROGRAM_H
