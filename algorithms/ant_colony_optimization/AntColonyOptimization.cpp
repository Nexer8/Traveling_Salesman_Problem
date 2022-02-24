//
// Created by Mariusz Wiśniewski on 2/23/2022.
//

#include <chrono>
#include <climits>
#include <cmath>
#include "AntColonyOptimization.hpp"

void AntColonyOptimization::calculateAntRoutes(Ant *ant, vector<vector<int>> &routes,
                                               vector<vector<double>> &pheromones) {
    vector<double> probabilities;

    routes[ant->idx][0] = ant->idx; // each ant's starting city depends on its index
    ant->visited[ant->idx] = true; // we mark the starting city as visited

    for (int i = 0; i < dimension - 1; ++i) {
        int city = routes[ant->idx][i]; // city of index i
        probabilities.clear();
        probabilities.resize(dimension, 0.0);
        for (int second_city = 0;
             second_city < dimension; second_city++) { // calculating the probabilities for all the cities
            if (city == second_city)
                continue;
            if (!ant->visited[second_city]) { // if the ant hasn't visited the city yet, then we calculate the probability
                probabilities[second_city] = calculateProbability(city, second_city, ant, pheromones);
            }
        }
        routes[ant->idx][i + 1] = getNextCity(probabilities); // deciding where to go next
        ant->visited[routes[ant->idx][i + 1]] = true; // marking the city as visited
    }
}

void AntColonyOptimization::dasAntColonyOptimization() {
    int result = INT_MAX;

    int number_of_ants = dimension; // we take the same amount of ants as the amount of the cities
    vector<vector<int>> ant_routes(number_of_ants);
    vector<vector<double>> pheromones(dimension);

    chrono::steady_clock::time_point begin = chrono::steady_clock::now();

    for (int i = 0; i < number_of_ants; ++i) {
        ant_routes[i].resize(number_of_ants, -1);
    }

    for (int i = 0; i < dimension; i++) {
        pheromones[i].resize(dimension);
        for (int j = 0; j < dimension; j++) {
            // setting up the same amount of pheromones at the beginning in order to avoid taking the greedy path
            pheromones[i][j] = (double) m_rand(0, RAND_MAX) / RAND_MAX * dimension / matrix[0][1];
        }
    }

    for (int i = 0; i < ACO_NUMBER_OF_ITERATIONS; i++) {
        for (int j = 0; j < number_of_ants; j++) {
            for (int &it: ant_routes[j]) {
                it = -1; // preparing the route
            }
            Ant *ant = new Ant(j, dimension); // creating a new ant
            calculateAntRoutes(ant, ant_routes, pheromones);
        }
        DAS(pheromones, ant_routes); // updating the pheromones values
    }


    for (int i = 0; i < dimension; i++) {
        int temp_res = calculateCost(ant_routes[i]);
        if (temp_res < result) {
            result = temp_res;
            best_path = ant_routes[i]; // choosing the best result and the best path
        }
    }

    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    cout << "Best DAS Ant Colony Optimization result: " << result << "\n";
    printSummary(begin, end);
}

void AntColonyOptimization::qasAntColonyOptimization() {
    int result = INT_MAX;

    int number_of_ants = dimension; // we take the same amount of ants as the amount of the cities
    vector<vector<int>> ant_routes(number_of_ants);
    vector<vector<double>> pheromones(dimension);

    chrono::steady_clock::time_point begin = chrono::steady_clock::now();

    for (int i = 0; i < number_of_ants; ++i) {
        ant_routes[i].resize(number_of_ants, -1);
    }

    for (int i = 0; i < dimension; i++) {
        pheromones[i].resize(dimension);
        for (int j = 0; j < dimension; j++) {
            // setting up the same amount of pheromones at the beginning in order to avoid taking the greedy path
            pheromones[i][j] = (double) m_rand(0, RAND_MAX) / RAND_MAX * dimension / matrix[0][1];
        }
    }

    for (int i = 0; i < ACO_NUMBER_OF_ITERATIONS; i++) {
        for (int j = 0; j < number_of_ants; j++) {
            for (int &it: ant_routes[j]) {
                it = -1; // preparing the route
            }
            Ant *ant = new Ant(j, dimension); // creating a new ant
            calculateAntRoutes(ant, ant_routes, pheromones);
        }
        QAS(pheromones, ant_routes); // updating the pheromones values
    }


    for (int i = 0; i < dimension; i++) {
        int temp_res = calculateCost(ant_routes[i]);
        if (temp_res < result) {
            result = temp_res;
            best_path = ant_routes[i]; // choosing the best result and the best path
        }
    }

    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    cout << "Best QAS Ant Colony Optimization result: " << result << "\n";
    printSummary(begin, end);
}

void AntColonyOptimization::casAntColonyOptimization() {
    int result = INT_MAX;

    int number_of_ants = dimension; // we take the same amount of ants as the amount of the cities
    vector<vector<int>> ant_routes(number_of_ants);
    vector<vector<double>> pheromones(dimension);

    chrono::steady_clock::time_point begin = chrono::steady_clock::now();

    for (int i = 0; i < number_of_ants; ++i) {
        ant_routes[i].resize(number_of_ants, -1);
    }

    for (int i = 0; i < dimension; i++) {
        pheromones[i].resize(dimension);
        for (int j = 0; j < dimension; j++) {
            // setting up the same amount of pheromones at the beginning in order to avoid taking the greedy path
            pheromones[i][j] = (double) m_rand(0, RAND_MAX) / RAND_MAX * dimension / matrix[0][1];
        }
    }

    for (int i = 0; i < ACO_NUMBER_OF_ITERATIONS; i++) {
        for (int j = 0; j < number_of_ants; j++) {
            for (int &it: ant_routes[j]) {
                it = -1; // preparing the route
            }
            Ant *ant = new Ant(j, dimension); // creating a new ant
            calculateAntRoutes(ant, ant_routes, pheromones);
        }
        CAS(pheromones, ant_routes); // updating the pheromones values
    }


    for (int i = 0; i < dimension; i++) {
        int temp_res = calculateCost(ant_routes[i]);
        if (temp_res < result) {
            result = temp_res;
            best_path = ant_routes[i]; // choosing the best result and the best path
        }
    }

    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    cout << "Best CAS Ant Colony Optimization result: " << result << "\n";
    printSummary(begin, end);
}

void AntColonyOptimization::CAS(vector<vector<double>> &pheromones, vector<vector<int>> &routes) {
    double q = dimension; // pheromone amount left on the route
    double ro = 0.5; // determines the amount of pheromone that evaporates in each iteration

    for (int i = 0; i < routes.size(); i++) { // for each route repeat
        int route = calculateCost(routes[i]); // calculate ant i's route cost
        for (int j = 0; j < routes.size() - 1; j++) {
            int city = routes[i][j]; // j city from i ant
            int next_city = routes[i][j + 1];

            // updating pheromone values on the edges between two cities
            pheromones[city][next_city] = (1 - ro) * pheromones[city][next_city] + q /
                                                                                   (double) route; // pheromone's amount decreases depending on the cost of a route
            pheromones[next_city][city] =
                    (1 - ro) * pheromones[next_city][city] + q / (double) route; // works in both ways
        }
    }
}

void AntColonyOptimization::QAS(vector<vector<double>> &pheromones, vector<vector<int>> &routes) {
    double q = dimension; // pheromone amount left on the route
    double ro = 0.5; // determines the amount of pheromone that evaporates in each iteration

    for (int i = 0; i < routes.size(); i++) { // for each route repeat
        for (int j = 0; j < routes.size() - 1; j++) {
            int city = routes[i][j]; // j city from i ant
            int next_city = routes[i][j + 1];

            // updating pheromone values on the edges between two cities
            pheromones[city][next_city] = (1 - ro) * pheromones[city][next_city] + q /
                                                                                   (double) matrix[city][next_city]; // pheromone's amount decreases depending on the cost of the edge
            pheromones[next_city][city] =
                    (1 - ro) * pheromones[next_city][city] + q / (double) matrix[next_city][city]; // works in both ways
        }
    }
}

void AntColonyOptimization::DAS(vector<vector<double>> &pheromones, vector<vector<int>> &routes) const {
    double q = dimension; // pheromone amount left on the route
    double ro = 0.5; // determines the amount of pheromone that evaporates in each iteration

    for (int i = 0; i < routes.size(); i++) { // for each route repeat
        for (int j = 0; j < routes.size() - 1; j++) {
            int city = routes[i][j]; // j city from i ant
            int next_city = routes[i][j + 1];

            // updating pheromone values on the edges between two cities
            pheromones[city][next_city] =
                    (1 - ro) * pheromones[city][next_city] + q; // pheromone's amount decreases depending on a constant
            pheromones[next_city][city] = (1 - ro) * pheromones[next_city][city] + q; // works in both ways
        }
    }
}

double AntColonyOptimization::calculateProbability(int first_city, int second_city, Ant *ant,
                                                   vector<vector<double>> &pheromones) {
    double alpha = 1; // alpha parameter adjusting the influence of pheromones
    double beta = 4.5; // beta parameter adjusting the influence of visibility

    // eta of going to another city
    auto eta_ij = (double) pow(1.0 / matrix[first_city][second_city], beta);
    auto tau_ij = (double) pow(pheromones[first_city][second_city], alpha);
    double sum = 0;

    for (int i = 0; i < dimension; ++i) {
        if (i == first_city) continue;

        if (!ant->visited[i]) {
            auto eta = (double) pow(1.0 / matrix[first_city][i], beta);
            auto tau = (double) pow(pheromones[first_city][i], alpha);
            sum += eta * tau; // a sum of what's left from the available cities
        }
    }

    return (eta_ij * tau_ij) / (sum);
}

void AntColonyOptimization::run() {
    switch (type) {
        case AntColonyOptimizationType::DENSITY:
            dasAntColonyOptimization();
            break;
        case AntColonyOptimizationType::QUANTITY:
            qasAntColonyOptimization();
            break;
        case AntColonyOptimizationType::CYCLE:
            casAntColonyOptimization();
            break;
    }
}

AntColonyOptimization::AntColonyOptimization(const Algorithm &alg, AntColonyOptimizationType type)
        : Algorithm(alg), type(type) {}
