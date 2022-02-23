//
// Created by Mariusz Wiśniewski on 2/23/2022.
//

#include <climits>
#include <chrono>
#include "TreeSearchBruteForce.hpp"

void TreeSearchBruteForce::treeSearch(vector<int> vertex_list) {
    if (vertex_list.size() == dimension) { // when the path is complete we check its cost
        vector<int> vec(dimension);
        int k = 0;
        for (int const &i: vertex_list) {
            vec[k++] = i;
        }
        int temp_result = calculateCost(vec);
        if (temp_result < min_cost) {
            min_cost = temp_result;
            for (int i = 0; i < dimension; i++) {
                best_path[i] = vec[i];
            }
            return;
        }
    }

    for (int i = 0; i < dimension; i++) { // checking which vertex we can add
        bool exists = false;
        for (int j: vertex_list) {
            if (i == j) {
                exists = true;
                break;
            }
        }
        if (!exists) { // if we can add a vertex we do this and continue the recursion
            vertex_list.push_back(i);
            treeSearch(vertex_list);
            vertex_list.pop_back(); // then we take the vertex out and check a new path
        }
    }
}

void TreeSearchBruteForce::run() {
    best_path = vector<int>(getDimension());
    chrono::steady_clock::time_point begin = chrono::steady_clock::now();

    vector<int> vertex_vector;
    min_cost = INT_MAX;
    vertex_vector.push_back(0);
    treeSearch(vertex_vector);
    cout << "Minimal cost: " << min_cost << "\n";

    chrono::steady_clock::time_point end = chrono::steady_clock::now();

    printSummary(begin, end);
}

TreeSearchBruteForce::TreeSearchBruteForce(const Algorithm &alg) : Algorithm(alg) {}
