//
// Created by mariu on 20/01/2020.
//
#include <vector>

#ifndef PEA_0_ANT_H
#define PEA_0_ANT_H

class Ant {
public:
    int idx;
    std::vector<bool> visited;
    Ant(int index, int number_of_vertex);
    ~Ant();
};


#endif //PEA_0_ANT_H
