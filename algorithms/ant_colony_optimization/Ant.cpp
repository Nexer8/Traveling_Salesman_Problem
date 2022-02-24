//
// Created by Mariusz Wiśniewski on 20/01/2020.
//

#include "Ant.hpp"

Ant::Ant(int index, int number_of_vertex) {
    idx = index;
    visited.resize(number_of_vertex, false);
}

Ant::~Ant() {
    visited.clear();
}