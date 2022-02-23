//
// Created by Mariusz Wiśniewski on 2/23/2022.
//

#include <climits>
#include <chrono>
#include "SwapBruteForce.hpp"

void SwapBruteForce::run() {
    vector<int> vec(getDimension());
    best_path = vector<int>(getDimension());
    for (int i = 0; i < getDimension(); i++) {
        vec[i] = i;
    }
    min_cost = INT_MAX;
    chrono::steady_clock::time_point begin_bf = chrono::steady_clock::now();
    swapBruteForce(0, vec);
    chrono::steady_clock::time_point end_bf = chrono::steady_clock::now();

    cout << "Minimal cost: " << min_cost << "\n";
    printSummary(begin_bf, end_bf);
}

void SwapBruteForce::swapBruteForce(int current_level, vector<int> &vec) {
    if (matrix.empty() || current_level < 0 || current_level > dimension) return;
    if (current_level == dimension) {
        int cost = calculateCost(vec);
        if (cost < min_cost) {
            min_cost = cost;
            for (int i = 0; i < dimension; i++) {
                best_path[i] = vec[i];
            }
        }
    }

    for (int i = current_level; i < dimension; i++) {
        // We swap the elements to generate a arr
        swap(vec, current_level, i);
        swapBruteForce(current_level + 1, vec);
        // Then we swap them back - backtracking
        swap(vec, current_level, i);
    }
}

SwapBruteForce::SwapBruteForce(const Algorithm &alg) : Algorithm(alg) {}
