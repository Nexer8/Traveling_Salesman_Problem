//
// Created by mariu on 20/01/2020.
//

#include "Ant.h"

Ant::Ant(int number, int number_of_vertex) {
    this->number = number;
    this->visited.resize(number_of_vertex, false);
}

Ant::~Ant() {
    this->visited.clear();
}