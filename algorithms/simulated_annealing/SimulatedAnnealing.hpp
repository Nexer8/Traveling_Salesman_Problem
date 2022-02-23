//
// Created by Mariusz Wiśniewski on 2/23/2022.
//

#ifndef TSP_SIMULATEDANNEALING_HPP
#define TSP_SIMULATEDANNEALING_HPP

#include "../../utils/constants.hpp"
#include "../../Algorithm.hpp"

#define STARTING_TEMPERATURE 1e9
#define STOPPING_TEMPERATURE 0.01
#define LINEAR_COOLING_CONST 0.95
#define GEOMETRIC_COOLING_CONST 5
#define REPEAT_VAL 3

#define NUMBER_OF_ITERATIONS 500

typedef double (*cooling_method)(double, int);

class SimulatedAnnealing : Algorithm {
private:
    static void chooseParams(vector<int> &vec, NeighborhoodType &nt, BeginningSolution &result, move_foo &move,
                             cooling_method &cooling, move_linear_update &linear_update);

public:
    void run() override;

    SimulatedAnnealing(Algorithm const & alg);
};


#endif //TSP_SIMULATEDANNEALING_HPP
