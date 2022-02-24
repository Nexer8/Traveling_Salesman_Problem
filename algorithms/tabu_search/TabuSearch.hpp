//
// Created by Mariusz Wiśniewski on 2/23/2022.
//

#ifndef TSP_TABUSEARCH_HPP
#define TSP_TABUSEARCH_HPP

#include "../Algorithm.hpp"

#define NUMBER_OF_ITERATIONS 500

#define STEPS_IN_GEN dimension
#define CADENCE dimension

class TabuSearch : Algorithm {
public:
    void run() override;

    explicit TabuSearch(Algorithm const &alg);
};


#endif //TSP_TABUSEARCH_HPP
