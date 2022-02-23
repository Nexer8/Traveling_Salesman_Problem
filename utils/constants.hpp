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

const string ATSP_PATH = "../data/ATSP/";
const string SMALL_PATH = "../data/SMALL/";
const string TSP_PATH = "../data/TSP/";

#endif //PEA_0_MAIN_PROGRAM_H
