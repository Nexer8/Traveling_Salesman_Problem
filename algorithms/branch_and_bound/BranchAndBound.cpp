//
// Created by Mariusz Wiśniewski on 2/23/2022.
//

#include <chrono>
#include <climits>
#include <cmath>
#include "BranchAndBound.hpp"

void BranchAndBound::run() {
    vector<int> tab = vector<int>(dimension);
    for (int i = 0; i < dimension; i++) tab[i] = i;

    int lower_bound = 0;
    // setting the lower bound so that is the sum of minimal distances for each vertex
    // We divide it by two because if not then we go through each path twice
    for (int i = 0; i < dimension; i++) {
        lower_bound += ceil((minDistanceFrom(i) + secondDistanceFrom(i)) / 2.0);
    }

    int level = 1;
    int distance = 0;
    min_cost = INT_MAX;
    best_path = vector<int>(dimension);
    chrono::steady_clock::time_point begin = chrono::steady_clock::now();

    branchNBound(tab, level, lower_bound, distance, min_cost);

    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    cout << "\nMinimal cost: " << min_cost << "\n";
    printSummary(begin, end);
}

void BranchAndBound::branchNBound(vector<int> &tab, int level, int lower_bound, int distance, int &min_distance) {
    if (level == dimension) {
        if (matrix[tab[level - 1]][tab[0]] < 0) return;

        // checking if the distance is smaller at the end
        distance += matrix[tab[level - 1]][tab[0]];
        if (distance < min_distance) {
            for (int i = 0; i < dimension; i++) best_path[i] = tab[i];
            min_distance = distance; // upper bound
        }
    } else {
        for (int i = level; i < dimension; i++) {
            if (matrix[tab[level - 1]][tab[i]] < 0) continue;

            // adjusting the lower_bound so that it fits the level we're currently at
            int temp_bound = lower_bound;
            if (level == 1) lower_bound -= ceil((minDistanceFrom(tab[0]) + minDistanceFrom(tab[i])) / 2.0);
            else if (level == dimension - 1)
                lower_bound -= ceil((secondDistanceFrom(tab[level - 1]) + secondDistanceFrom(tab[0])) / 2.0);
            else lower_bound -= ceil((secondDistanceFrom(tab[level - 1]) + minDistanceFrom(tab[i])) / 2.0);

            // distance from 0 to i
            distance += matrix[tab[level - 1]][tab[i]];

            // verify if further checking makes sense
            // lower_bound + actual distance should be then smaller than the minimal distance found, if not, we cut the node
            if (lower_bound + distance < min_distance) {
                swap(tab, level, i);

                branchNBound(tab, level + 1, lower_bound, distance, min_distance);

                swap(tab, level, i);
            }

            distance -= matrix[tab[level - 1]][tab[i]]; // we go back in a tree
            lower_bound = temp_bound;
        }
    }
}

int BranchAndBound::minDistanceFrom(int vertex) {
    int min = INT_MAX;
    for (int i = 0; i < dimension; i++) {
        if (matrix[vertex][i] < min && i != vertex)
            min = matrix[vertex][i]; // finding the minimum edge cost having an end at the vertex "vertex"
    }
    return min;
}

int BranchAndBound::secondDistanceFrom(int vertex) {
    int min = INT_MAX, second = INT_MAX;
    for (int i = 0; i < dimension; i++) {
        if (i == vertex) continue;

        if (matrix[vertex][i] < min && i != vertex) {
            second = min;
            min = matrix[vertex][i]; // finding the second minimum edge cost having an end at the vertex "vertex"
        } else if (matrix[vertex][i] <= second && matrix[vertex][i] != min) {
            second = matrix[vertex][i];
        }
    }
    return second;
}

BranchAndBound::BranchAndBound(const Algorithm &alg) : Algorithm(alg) {}
