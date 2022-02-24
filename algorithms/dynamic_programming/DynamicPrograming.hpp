//
// Created by Mariusz Wiśniewski on 2/23/2022.
//

#ifndef TSP_DYNAMICPROGRAMING_HPP
#define TSP_DYNAMICPROGRAMING_HPP

#include "../Algorithm.hpp"

class DynamicProgramming : Algorithm {
public:
    void run() override;

    explicit DynamicProgramming(Algorithm const &alg);
};


#endif //TSP_DYNAMICPROGRAMING_HPP
