//
// Created by Mariusz Wiśniewski on 2/23/2022.
//

#include <vector>
#include <chrono>
#include <climits>
#include "TabuSearch.hpp"
#include "../../utils/constants.hpp"


void TabuSearch::run() {
    move_foo move;
    move_linear_update linear_update;
    NeighborhoodType ngbh_type;

    best_path = vector<int>(dimension);
    BeginningSolution beg_solution;
    int result;
    vector<int> vec(dimension);
    for (int i = 0; i < dimension; i++) {
        vec[i] = i;
    }

    // Setting a type of the neighbourhood and the beginning solution
    chooseCommonParams(vec, ngbh_type, beg_solution, move, linear_update);
    chrono::steady_clock::time_point begin = chrono::steady_clock::now();

    switch (beg_solution) {
        case NATURAL:
            result = calculateCost(vec);
            break;
        case RANDOM:
            result = genRandomPathCost(vec);
            break;
        case GREEDY:
            result = greedy(vec, 0);
            break;
    }
    int current_val = result;
    vector<int> temp_path(dimension);

    // Creating a tabu_list
    vector<vector<int>> tabu_list;
    tabu_list.resize(dimension);
    for (auto &it: tabu_list) {
        it.resize(dimension, 0);
    }

    vector<int> next_step(dimension);

    for (int num_of_iteration = 0;
         num_of_iteration < NUMBER_OF_ITERATIONS; num_of_iteration++) { // setting a number of iterations
        for (int step = 0; step < STEPS_IN_GEN * dimension; step++) {   // number of steps in a generation
            int next_step_val = INT_MAX;

            int f_tabu = 0, s_tabu = 0;

            for (int i = 0; i < dimension; i++) {
                for (int j = i + 1; j < dimension; j++) {
                    int current_val_tmp = current_val;

                    if (ngbh_type == INSERT) {
                        temp_path = vec;
                    }
                    if ((j - i < 3) || (j == dimension - 1) ||
                        (i == 0) && ngbh_type != INSERT) { // boundary cases, just update in O(n) time
                        move(vec, i, j);
                        current_val = calculateCost(vec);
                    } else if ((i == 0 || j == dimension - 1) && ngbh_type == INSERT) {
                        move(vec, i, j);
                        current_val = calculateCost(vec);
                    } else {  // O(1) update
                        current_val = linear_update(current_val, vec, i, j, matrix);
                        move(vec, i, j);
                    }

                    // if the solution in neighbourhood is better than the current global best
                    if (current_val < result) {
                        result = current_val; // aspiration criteria - don't even check the tabu_list
                        best_path = vec;// copy solution array to best_path array
                    }

                    if (current_val <
                        next_step_val) {  // if the solution is better than the local best (inside a given step)
                        if (tabu_list[j][i] <=
                            step) {  // if the move is available, then update tabu_list and next_step_val
                            s_tabu = i;
                            f_tabu = j;
                            next_step = vec;
                            next_step_val = current_val;
                        }
                    }
                    if (ngbh_type != INSERT) move(vec, i, j); // go back to the previous position
                    else {
                        vec = temp_path;
                    }
                    current_val = current_val_tmp;
                }
            }
            vec = next_step;
            tabu_list[f_tabu][s_tabu] += CADENCE; // updating tabu list with the opposite value
            // tabu on a given edge exists for <dimension> iterations
        }

        // generate a new set - critical event
        switch (beg_solution) {
            case NATURAL:
                for (int i = 0; i < dimension; i++) {
                    vec[i] = i;
                }
                current_val = calculateCost(vec);
                break;
            case RANDOM:
                current_val = genRandomPathCost(vec);
                break;
            case GREEDY:
                current_val = greedy(vec, 0);
                break;
        }
        // restoring tabu_list
        for (auto &it: tabu_list) {
            it.clear();
        }
        tabu_list.clear();
        tabu_list.resize(dimension);
        for (auto &it: tabu_list) {
            it.resize(dimension, 0);
        }
    }
    chrono::steady_clock::time_point end = chrono::steady_clock::now();

    cout << "Best Tabu Search result: " << result << "\n";
    printSummary(begin, end);
}

TabuSearch::TabuSearch(const Algorithm &alg) : Algorithm(alg) {}
