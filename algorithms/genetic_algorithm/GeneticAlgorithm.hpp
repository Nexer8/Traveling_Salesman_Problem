//
// Created by Mariusz Wiśniewski on 2/23/2022.
//

#ifndef TSP_GENETICALGORITHM_HPP
#define TSP_GENETICALGORITHM_HPP

#include "../Algorithm.hpp"

#define POPULATION_SIZE 5000
#define NUMBER_OF_GENERATIONS 500
#define NUMBER_OF_TOURNAMENTS 5
#define CROSS_RATE 0.2

#define MUTATION_RATE 0.02

class GeneticAlgorithm : Algorithm {
private:
    static void chooseParams(move_foo &move, NeighborhoodType &nt, crossover_type &crossover);

    void generatePopulation();

    void selectMatingPoolTournament();

    static void OX(vector<int> &first_parent, vector<int> &second_parent);

    static void PMX(vector<int> &first_parent, vector<int> &second_parent);

    static void NWOX(vector<int> &first_parent, vector<int> &second_parent);

public:
    void run() override;

    explicit GeneticAlgorithm(Algorithm const &alg);
};


#endif //TSP_GENETICALGORITHM_HPP
