//
// Created by Mariusz Wiśniewski on 2/23/2022.
//

#ifndef TSP_ANTCOLONYOPTIMIZATION_HPP
#define TSP_ANTCOLONYOPTIMIZATION_HPP

#include "../../Algorithm.hpp"
#include "Ant.hpp"

#define ACO_NUMBER_OF_ITERATIONS 100

enum class AntColonyOptimizationType {
    DENSITY,
    QUANTITY,
    CYCLE
};

class AntColonyOptimization : Algorithm {
private:
    AntColonyOptimizationType type = AntColonyOptimizationType::DENSITY;

    void calculate_ant_routes(Ant *ant, vector<vector<int>> &routes, vector<vector<double>> &pheromones);

    void das_ant_colony_optimization();

    void cas_ant_colony_optimization();

    void qas_ant_colony_optimization();

    void CAS(vector<vector<double>> &pheromones, vector<vector<int>> &routes);

    void DAS(vector<vector<double>> &pheromones, vector<vector<int>> &routes) const;

    void QAS(vector<vector<double>> &pheromones, vector<vector<int>> &routes);

    double calculate_probability(int first_city, int second_city, Ant *ant, vector<vector<double>> &pheromones);

public:
    void run() override;

    AntColonyOptimization(Algorithm const &alg, AntColonyOptimizationType type);
};


#endif //TSP_ANTCOLONYOPTIMIZATION_HPP
