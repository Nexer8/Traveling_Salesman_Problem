//
// Created by Mariusz Wiśniewski on 2/23/2022.
//

#ifndef TSP_TREESEARCHBRUTEFORCE_HPP
#define TSP_TREESEARCHBRUTEFORCE_HPP

#include "../Algorithm.hpp"

class TreeSearchBruteForce : Algorithm {
private:
    void treeSearch(vector<int> vertex_list);

public:
    explicit TreeSearchBruteForce(Algorithm const &alg);

    void run() override;
};


#endif //TSP_TREESEARCHBRUTEFORCE_HPP
