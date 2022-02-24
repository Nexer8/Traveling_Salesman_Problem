//
// Created by Mariusz Wiśniewski on 2/23/2022.
//

#ifndef TSP_BRANCHANDBOUND_HPP
#define TSP_BRANCHANDBOUND_HPP


#include "../Algorithm.hpp"

class BranchAndBound : Algorithm {
private:
    void branchNBound(vector<int> &tab, int level, int lower_bound, int distance, int &min_distance);

    int minDistanceFrom(int v);

    int secondDistanceFrom(int v);

public:
    void run() override;

    explicit BranchAndBound(Algorithm const &alg);
};


#endif //TSP_BRANCHANDBOUND_HPP
