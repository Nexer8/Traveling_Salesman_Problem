//
// Created by Mariusz Wiśniewski on 2/23/2022.
//

#ifndef TSP_SWAPBRUTEFORCE_HPP
#define TSP_SWAPBRUTEFORCE_HPP

#include "../../Algorithm.hpp"

class SwapBruteForce : Algorithm {
private:
    void swapBruteForce(int current_level, vector<int> &vec);

public:
    explicit SwapBruteForce(Algorithm const & alg);

    void run() override;
};


#endif //TSP_SWAPBRUTEFORCE_HPP
