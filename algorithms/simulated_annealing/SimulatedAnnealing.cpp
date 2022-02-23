//
// Created by Mariusz Wiśniewski on 2/23/2022.
//

#include <cmath>
#include <chrono>
#include <cstring>
#include "SimulatedAnnealing.hpp"

static inline double linear_cooling(double temperature, int it) {
    return temperature *= LINEAR_COOLING_CONST;
}

static inline double geometric_cooling(double temperature, int it) {
    return temperature -= GEOMETRIC_COOLING_CONST;
}

static inline double logarithmic_cooling(double temperature, int it) {
    return STARTING_TEMPERATURE / (log(it) + 1);
}

static inline double get_probability(int difference, double temperature) { // checking how weak the solution is
    return exp(-1 * difference / temperature);
}

void SimulatedAnnealing::chooseParams(vector<int> &vec, NeighborhoodType &nt, BeginningSolution &result, move_foo &move,
                                      cooling_method &cooling, move_linear_update &linear_update) {
    chooseCommonParams(vec, nt, result, move, linear_update);
    char option;

    cout << "\n";
    cout << "Choose cooling method (linear by default): \n";
    cout << " 1 - Linear.\n";
    cout << " 2 - Geometric.\n";
    cout << " 3 - Logarithmic.\n";
    cout << "\n";
    cin >> option;

    switch (option) {
        case '1':
            cooling = linear_cooling;
            break;
        case '2':
            cooling = geometric_cooling;
            break;
        case '3':
            cooling = logarithmic_cooling;
            break;
        default:
            cooling = linear_cooling;
            break;
    }
}

void SimulatedAnnealing::run() {
    move_foo move;
    cooling_method cooling;
    move_linear_update linear_update;
    NeighborhoodType ngbh_type;

    best_path = vector<int>(dimension);
    vector<int> vec(dimension);
    for (int i = 0; i < dimension; i++) {
        vec[i] = i;
    }
    BeginningSolution beg_solution;
    int result;

    chooseParams(vec, ngbh_type, beg_solution, move, cooling, linear_update);

    chrono::steady_clock::time_point begin = chrono::steady_clock::now();

    switch (beg_solution) {
        case NATURAL:
            result = calculateCost(vec);
            break;
        case RANDOM:
            result = genRandomPathCost(vec);
            break;
        case GREEDY:
            result = greedy(vec, m_rand(0, dimension - 1));
            break;
    }

    int step_val = result;

    vector<int> next_step(dimension);
    vector<int> temp_arr(dimension);

    int stagnation = 0; // stagnation value
    int iteration = 0;  // number of iterations for one annealing and cooling
    double temperature = STARTING_TEMPERATURE;

    for (int num_of_iteration = 0; num_of_iteration < NUMBER_OF_ITERATIONS; ++num_of_iteration) {
        while (temperature > STOPPING_TEMPERATURE) { // end condition
            iteration++;
            int repeat = REPEAT_VAL *
                         dimension;  // number of neighbourhood verifications = 3n, number of iterations for the same temperature
            int temp_step_val = step_val;

            while (repeat-- > 0) {
                int i = 0;  // first iterator
                int j = 0;  // second iterator
                while (i == j) {
                    i = m_rand(0, dimension - 1);
                    j = m_rand(0, dimension - 1);
                }
                if (i > j) ::swap(i, j);

                if (ngbh_type == INSERT) {
                    temp_arr = vec;
                }
                // checking step_value, O(n) update
                if ((j - i < 3) || j == dimension - 1 || i == 0 && ngbh_type != INSERT) { // border cases verification
                    move(vec, i, j);
                    step_val = calculateCost(vec);
                } else if ((i == 0 || j == dimension - 1) && ngbh_type == INSERT) {
                    move(vec, i, j);
                    step_val = calculateCost(vec);
                } // O(1) update
                else {
                    step_val = linear_update(step_val, vec, i, j, matrix);
                    move(vec, i, j);
                }

                int difference = step_val - result; // calculate difference between the best and actual results

                if (step_val <
                    result) {  // if step_val is better than current result, update result and reset stagnation value
                    stagnation = 0;
                    result = step_val;
                    best_path = vec;
                } else {   // if the value has not improved, increase the stagnation value
                    stagnation++;
                    if (stagnation > dimension)
                        temperature = (double) (STARTING_TEMPERATURE) /
                                      2;   // if stagnation is fairly big, increase the temperature
                }

                // decision-making whether to take the next_step
                if (difference < 0 ||
                    (difference > 0 &&
                     get_probability(difference, temperature) > ((double) m_rand(0, RAND_MAX) / RAND_MAX))) {
                    next_step = vec;
                    break;
                } else {
                    if (ngbh_type != INSERT) move(vec, i, j); // go back to the previous position
                    else {
                        vec = temp_arr;
                    }
                    step_val = temp_step_val;
                }
            }
            temperature = cooling(temperature, iteration);  // cooling
            vec = next_step; // assigning best step
            if (iteration > dimension * dimension) break; // if we are for too long in the same spot
        }

        // generating new set of starting conditions
        switch (beg_solution) {
            case NATURAL:
                for (int i = 0; i < dimension; i++) {
                    vec[i] = i;
                }
                step_val = calculateCost(vec);
                break;
            case RANDOM:
                step_val = genRandomPathCost(vec);
                break;
            case GREEDY:
                step_val = greedy(vec, m_rand(0, dimension - 1));
                break;
        }
        iteration = 0;
        stagnation = 0;
        temperature = STARTING_TEMPERATURE;
    }
    chrono::steady_clock::time_point end = chrono::steady_clock::now();

    cout << "Best Simulated Annealing result: " << result << "\n";
    printSummary(begin, end);
}

SimulatedAnnealing::SimulatedAnnealing(const Algorithm &alg) : Algorithm(alg) {}
