//
// Created by Mariusz Wiśniewski on 2/23/2022.
//

#include <algorithm>
#include <chrono>
#include <climits>
#include "GeneticAlgorithm.hpp"

void GeneticAlgorithm::OX(vector<int> &first_parent, vector<int> &second_parent) {
    int dimension = static_cast<int>(first_parent.size());
    int k1, k2;
    vector<int> first_child(dimension, -1);
    vector<int> second_child(dimension, -1);

    do {
        k1 = m_rand(0, dimension - 2);
        k2 = m_rand(0, dimension - 2);  //in case of hitting last index, while loop won't break (numbers from 1 to n-1)
    } while (k1 == k2);

    if (k1 > k2) {
        ::swap(k1, k2);
    }

    for (int i = k1; i <= k2; i++) {
        first_child[i] = second_parent[i];
        second_child[i] = first_parent[i];
    }

    auto child_iterator = first_child.begin() + k2 + 1;
    auto parent_iterator = first_parent.begin() + k2 + 1;

    while (child_iterator != first_child.begin() + k1) {
        if (first_child.end() ==
            find(first_child.begin(), first_child.end(),
                 *parent_iterator)) { // if a child does not have parent's gene (city)
            *child_iterator = *parent_iterator; // value of child's gene = value of parent's gene

            if (child_iterator == first_child.end() - 1)    // when child's iterator achieves the last element
                child_iterator = first_child.begin();       // go to the beginning
            else
                child_iterator++;

            if (parent_iterator == first_parent.end() - 1) // same for the parent
                parent_iterator = first_parent.begin();
            else
                parent_iterator++;
        } else { // if a child already has that gene (city), go to the next parent's gene (city)
            if (parent_iterator == first_parent.end() - 1)
                parent_iterator = first_parent.begin();
            else
                parent_iterator++;
        }
    }

    child_iterator = second_child.begin() + k2 + 1;
    parent_iterator = second_parent.begin() + k2 + 1;

    // second child follows suit
    while (child_iterator != second_child.begin() + k1) {
        if (second_child.end() == find(second_child.begin(), second_child.end(), *parent_iterator)) {
            *child_iterator = *parent_iterator;

            if (child_iterator == second_child.end() - 1)
                child_iterator = second_child.begin();
            else
                child_iterator++;

            if (parent_iterator == second_parent.end() - 1)
                parent_iterator = second_parent.begin();
            else
                parent_iterator++;
        } else {
            if (parent_iterator == second_parent.end() - 1)
                parent_iterator = second_parent.begin();
            else
                parent_iterator++;
        }
    }

    first_parent = first_child;
    second_parent = second_child;
}

void GeneticAlgorithm::PMX(vector<int> &first_parent, vector<int> &second_parent) {
    int dimension = static_cast<int>(first_parent.size());
    int k1, k2;
    vector<int> first_child(dimension, -1);
    vector<int> second_child(dimension, -1);

    struct pair {
        int first_idx;
        int second_idx;
    };

    vector<pair> map;

    do {
        k1 = m_rand(0, dimension - 2);
        k2 = m_rand(0, dimension - 2);  //in case of hitting last index, while loop won't break (numbers from 1 to n-1)
    } while (k1 == k2);

    if (k1 > k2) {
        ::swap(k1, k2);
    }

    for (int i = k1; i <= k2; i++) {
        first_child[i] = second_parent[i];
        second_child[i] = first_parent[i];
        pair p = {.first_idx = first_parent[i], .second_idx = second_parent[i]};
        map.push_back(p);
    }

    // First child
    for (int i = 0; i < dimension; i++) {
        if (first_child[i] != -1) continue;

        bool exists_in_map = false;
        for (auto &j: map) {
            if (first_parent[i] == j.first_idx || first_parent[i] == j.second_idx) {
                exists_in_map = true;
                break;
            }
        }

        if (first_child[i] == -1
            && first_child.end() == find(first_child.begin(), first_child.end(), first_parent[i])
            && !exists_in_map) {
            first_child[i] = first_parent[i];
        }
    }

    for (int i = 0; i < dimension; i++) {
        int val = first_parent[i];
        if (first_child[i] == -1) {
            vector<int> used_values;
            used_values.push_back(val);
            while (first_child.end() != find(first_child.begin(), first_child.end(), val)) {
                for (auto &j: map) {
                    if (val == j.first_idx &&
                        used_values.end() == find(used_values.begin(), used_values.end(), j.second_idx)) {
                        val = j.second_idx;
                        used_values.push_back(val);
                        break;
                    } else if (val == j.second_idx &&
                               used_values.end() == find(used_values.begin(), used_values.end(), j.first_idx)) {
                        val = j.first_idx;
                        used_values.push_back(val);
                        break;
                    }
                }
            }
            first_child[i] = val;
        }
    }

    // Second child
    for (int i = 0; i < dimension; i++) {
        if (second_child[i] != -1) continue;
        bool exists_in_map = false;
        for (auto &j: map) {
            if (second_parent[i] == j.first_idx || second_parent[i] == j.second_idx) {
                exists_in_map = true;
                break;
            }
        }

        if (second_child[i] == -1
            && second_child.end() == find(second_child.begin(), second_child.end(), second_parent[i])
            && !exists_in_map) {
            second_child[i] = second_parent[i];
        }
    }

    for (int i = 0; i < dimension; i++) {
        int val = second_parent[i];
        if (second_child[i] == -1) {
            vector<int> used_values;
            used_values.push_back(val);
            while (second_child.end() != find(second_child.begin(), second_child.end(), val)) {
                for (auto &j: map) {
                    if (val == j.first_idx &&
                        used_values.end() == find(used_values.begin(), used_values.end(), j.second_idx)) {
                        val = j.second_idx;
                        used_values.push_back(val);
                        break;
                    } else if (val == j.second_idx &&
                               used_values.end() == find(used_values.begin(), used_values.end(), j.first_idx)) {
                        val = j.first_idx;
                        used_values.push_back(val);
                        break;
                    }
                }
            }
            second_child[i] = val;
        }
    }
    first_parent = first_child;
    second_parent = second_child;
}

void GeneticAlgorithm::NWOX(vector<int> &first_parent, vector<int> &second_parent) {
    int dimension = static_cast<int>(first_parent.size());
    int k1, k2;
    vector<int> first_child(dimension, -1);
    vector<int> second_child(dimension, -1);

    do {
        k1 = m_rand(0, dimension - 2);
        k2 = m_rand(0, dimension - 2);  //in case of hitting last index, while loop won't break (numbers from 1 to n-1)
    } while (k1 == k2);

    if (k1 > k2) {
        ::swap(k1, k2);
    }

    for (int i = k1; i <= k2; i++) {
        first_child[i] = second_parent[i];
        second_child[i] = first_parent[i];
    }

    // First child
    int first_child_dimension = first_child.size();

    for (int i = 0; i < first_child_dimension; i++) {
        if (first_child[i] == -1 &&
            first_child.end() == find(first_child.begin(), first_child.end(), first_parent[i])) {
            first_child[i] = first_parent[i];
        }
    }

    for (int i = 0; i < k1 - 1; i++) {
        if (first_child[i] == -1) {
            for (int j = i + 1; j < k1; j++) {
                if (first_child[j] != -1) {
                    first_child[i] = first_child[j];
                    first_child[j] = -1;
                    break;
                }
            }
        }
    }

    for (int i = first_child_dimension - 1; i > k2; i--) {
        if (first_child[i] == -1) {
            for (int j = i - 1; j > k2; j--) {
                if (first_child[j] != -1) {
                    first_child[i] = first_child[j];
                    first_child[j] = -1;
                    break;
                }
            }
        }
    }

    for (int i = 0; i < first_child_dimension; i++) {
        if (first_child[i] == -1) {
            for (int j = k1; j <= k2; j++) {
                if (first_child[i] == -1 &&
                    first_child.end() == find(first_child.begin(), first_child.end(), first_parent[j])) {
                    first_child[i] = first_parent[j];
                    break;
                }
            }
        }
    }

    // Second child
    int second_child_dimension = second_child.size();

    for (int i = 0; i < second_child_dimension; i++) {
        if (second_child[i] == -1 &&
            second_child.end() == find(second_child.begin(), second_child.end(), second_parent[i])) {
            second_child[i] = second_parent[i];
        }
    }

    for (int i = 0; i < k1 - 1; i++) {
        if (second_child[i] == -1) {
            for (int j = i + 1; j < k1; j++) {
                if (second_child[j] != -1) {
                    second_child[i] = second_child[j];
                    second_child[j] = -1;
                    break;
                }
            }
        }
    }

    for (int i = second_child_dimension - 1; i > k2; i--) {
        if (second_child[i] == -1) {
            for (int j = i - 1; j > k2; j--) {
                if (second_child[j] != -1) {
                    second_child[i] = second_child[j];
                    second_child[j] = -1;
                    break;
                }
            }
        }
    }

    for (int i = 0; i < second_child_dimension; i++) {
        if (second_child[i] == -1) {
            for (int j = k1; j <= k2; j++) {
                if (second_child[i] == -1 &&
                    second_child.end() == find(second_child.begin(), second_child.end(), second_parent[j])) {
                    second_child[i] = second_parent[j];
                    break;
                }
            }
        }
    }

    first_parent = first_child;
    second_parent = second_child;
}


void GeneticAlgorithm::chooseParams(move_foo &move, NeighborhoodType &nt, crossover_type &crossover) {
    char option;

    cout << "\n";
    cout << "Choose mutation type (swap by default): \n";
    cout << " 1 - Swap.\n";
    cout << " 2 - Invert.\n";
    cout << " 3 - Insert.\n";
    cout << "\n";
    cin >> option;

    switch (option) {
        case '1':
            move = swap;
            nt = SWAP;
            break;
        case '2':
            move = invert;
            nt = INVERT;
            break;
        case '3':
            move = insert;
            nt = INSERT;
            break;
        default:
            move = swap;
            break;
    }

    cout << "\n";
    cout << "Choose crossover type (OX by default): \n";
    cout << " 1 - OX.\n";
    cout << " 2 - PMX.\n";
    cout << " 3 - NWOX.\n";
    cout << "\n";
    cin >> option;

    switch (option) {
        case '1':
            crossover = OX;
            break;
        case '2':
            crossover = PMX;
            break;
        case '3':
            crossover = NWOX;
            break;
        default:
            crossover = OX;
            break;
    }
}

void GeneticAlgorithm::run() {
    generatePopulation();
    int result = INT_MAX;

    move_foo mutate; // mutation
    NeighborhoodType ngbh_type;
    crossover_type crossover; // crossover operator - OX, PMX and NWOX available

    chooseParams(mutate, ngbh_type, crossover); // choose mutation and crossover types
    chrono::steady_clock::time_point begin = chrono::steady_clock::now();

    best_path = population[0]; // at first, the best path is the first available one
    for (int j = 0; j < NUMBER_OF_GENERATIONS; j++) { // repeat all for number of generations
        int first_rand_idx;
        int second_rand_idx;
        selectMatingPoolTournament(); // select the mating pool
        for (int i = 0;
             i < (int) (POPULATION_SIZE * CROSS_RATE); i++) { // CROSS_RATE tells us how often we do crossovers
            do {
                first_rand_idx = m_rand(0, POPULATION_SIZE - 1);
                second_rand_idx = m_rand(0, POPULATION_SIZE - 1);
            } while (first_rand_idx == second_rand_idx);

            crossover(population[first_rand_idx], population[second_rand_idx]); // the actual crossover operation
        }

        for (int i = 0; i < (int) (POPULATION_SIZE * MUTATION_RATE); i++) { // then we do the mutation
            int rand_index = m_rand(0, POPULATION_SIZE - 1);
            do {
                first_rand_idx = m_rand(0, dimension - 1);
                second_rand_idx = m_rand(0, dimension - 1);
            } while (first_rand_idx == second_rand_idx);
            mutate(population[rand_index], first_rand_idx, second_rand_idx); // the actual mutation operation
        }

        for (int i = 0; i <
                        POPULATION_SIZE; i++) { // choosing which population has the best fitness function - which path has the lowest cost
            fitness[i] = calculateCost(population[i]);
            if (result > fitness[i]) {
                result = fitness[i];
                best_path = population[i]; // getting the best path and the best solution
            }
        }
    }

    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    cout << "Best Genetic Algorithm result: " << result << "\n";
    printSummary(begin, end);
}

void GeneticAlgorithm::generatePopulation() {
    population.clear();
    fitness.clear();
    population.reserve(POPULATION_SIZE);
    fitness.reserve(POPULATION_SIZE);

    vector<int> route;
    for (int i = 0; i < dimension; i++) {
        route.push_back(i);
    }

    // draw new permutations and push them to the population vector
    for (int i = 0; i < POPULATION_SIZE; ++i) {
        fitness.push_back(genRandomPathCost(route));
        population.push_back(route);
    }
}

void GeneticAlgorithm::selectMatingPoolTournament() { // tournament mating pool selection
    vector<vector<int>> mating_pool;
    mating_pool.reserve(POPULATION_SIZE);

    for (int j = 0; j < POPULATION_SIZE; j++) {
        int best = INT_MAX;
        int best_idx = 0;
        for (int i = 0; i < NUMBER_OF_TOURNAMENTS; i++) {
            int rand_idx = m_rand(0, POPULATION_SIZE - 1); // draw a random index
            int rand_idx_val = fitness[rand_idx];
            if (best > rand_idx_val) { // choose the best permutation from the randomly generated ones
                best = rand_idx_val;
                best_idx = rand_idx;
            }
        }
        mating_pool.push_back(population[best_idx]);
    }
    population = mating_pool;
    for (int i = 0; i < POPULATION_SIZE; i++) {
        fitness[i] = calculateCost(population[i]);
    }
}

GeneticAlgorithm::GeneticAlgorithm(const Algorithm &alg) : Algorithm(alg) {}
