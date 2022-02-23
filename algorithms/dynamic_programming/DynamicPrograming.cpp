//
// Created by Mariusz Wiśniewski on 2/23/2022.
//

#include <climits>
#include <chrono>
#include "DynamicPrograming.hpp"

void DynamicProgramming::run() {
    best_path = vector<int>(getDimension());
    min_cost = INT_MAX;

    chrono::steady_clock::time_point begin_dp = chrono::steady_clock::now();

    vector<vector<int>> dp_matrix(1u << (unsigned) dimension); // [path][last vertex]
    // create a room for all possible paths

    for (unsigned int i = 0; i < (1u << (unsigned) dimension); i++) {
        dp_matrix[i] = vector<int>(dimension); // number of columns is equal to matrix dimension
    }

    for (unsigned int i = 0; i < (1u << (unsigned) dimension); i++) {
        for (int j = 0; j < dimension; j++) {
            dp_matrix[i][j] = 1u << 30u; // set all the costs to "infinity"
        }
    }

    for (unsigned int i = 0; i < dimension; i++) {
        if (i == 0) dp_matrix[1u << i | 1u][i] = 0; // Don't want a -1 value to be put in here.
        else {
            dp_matrix[1u << i | 1u][i] = matrix[0][i]; // a cost of 0<->(n-1) path
        }
    }

    for (unsigned int bit_mask = 0; bit_mask < 1u << (unsigned) dimension; bit_mask++) {
        for (unsigned int vertex = 0; vertex < (unsigned) dimension; vertex++) {
            if (!(bit_mask & (1u << vertex)))
                continue; // filter just the paths that don't include the vertex and has to end in it
            for (unsigned int j = 0; j < (unsigned) dimension; j++) {
                if (!(bit_mask &
                      (1u << j))) {// vertex j which is a proposed next extension cannot already exist in the path
                    if (dp_matrix[bit_mask][vertex] + matrix[vertex][j] < dp_matrix[bit_mask | (1u
                            << j)][j]) { // check if it's better to add vertex j to the current path that finishes in vertex or the actual path till j is more optimal
                        dp_matrix[bit_mask | (1u << j)][j] = dp_matrix[bit_mask][vertex] + matrix[vertex][j];
                    }
                }
            }
        }
    }

    int cost;
    for (int i = 0; i < dimension; i++) { // looking for the minimal cost cycle in the last row
        cost = dp_matrix[(1u << (unsigned) dimension) - 1][i] + matrix[i][0];
        if (cost < min_cost) {
            min_cost = cost;
        }
    }

    cout << "Minimal cost: " << min_cost << "\n";
    chrono::steady_clock::time_point end_dp = chrono::steady_clock::now();
    printDuration(begin_dp, end_dp);
}

DynamicProgramming::DynamicProgramming(const Algorithm &alg) : Algorithm(alg) {}
