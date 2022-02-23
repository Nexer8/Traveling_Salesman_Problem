//
// Created by Mariusz Wiśniewski on 20/01/2020.
//

#include "Ant.hpp"

Ant::Ant(int index, int number_of_vertex) {
    this->idx = index;
    this->visited.resize(number_of_vertex, false);
}

Ant::~Ant() {
    this->visited.clear();
}