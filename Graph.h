#include <iostream>
#include <string>
#include <vector>
#include <c++/4.8.3/functional>
#include "Ant.h"

#ifndef PEA_0_GRAPH_H
#define PEA_0_GRAPH_H

#define EXEC_ERROR -1

using namespace std;

typedef void (*move_foo)(int*, int, int);
typedef double (*cooling_method)(double, int);
typedef int (*move_linear_update)(int*, int*, int*, int*, int**);

enum beginning_solution {
    NATURAL,
    RANDOM,
    GREEDY
};

enum neighbourhood_type {
    SWAP,
    INVERT,
    INSERT
};

enum problem_type {
    TSP,
    SMALL,
    ATSP
};

class Graph {
private:
    int **matrix;       // the actual graph saved as a matrix
    int dimension;      // dimension of the graph

public:
    int *best_path;
    int min_cost;
    vector<vector< int > > population;
    vector<int> fitness;

    Graph();
    ~Graph();
    int load_data(problem_type version, string file_name);
    int print_graph();
    int calculate_cost(const int *arr);
    int get_dimension();
    void brute_force_swap(int current_level, int *arr);
    void dynamic_programming();
    void branch_n_bound(int *tab, int level, int bound, int distance, int& minDistance);
    void branch_n_bound();
    int min_distance_from(int v);
    int second_distance_from(int v);
    void tree_search(vector<int> vertex_list);
    void brute_force_tree();
    int greedy(int *arr, int starting_point);
    int gen_random_path_cost(int *arr);
    void tabu_search();
    void simulated_annealing();
    void ts_choose_params(int *arr, neighbourhood_type *nt, beginning_solution *result, move_foo *move, move_linear_update *linear_update);
    void sa_choose_params(int *arr, neighbourhood_type *nt, beginning_solution *result, move_foo *move, cooling_method *cooling, move_linear_update *linear_update);
    void generate_population();
    vector<int> gen_random_path();
    void choose_population();
    void ordered_crossover(vector<int> &first_parent, vector<int> &second_parent);
    void inversion_mutation(vector<int> &path);
    void update_pheromones(vector<vector<double>> &pheromones, vector<vector<int>> &routes);
    double get_phi(int first_city, int second_city, Ant *ant, vector<vector<double>> &pheromones);
    int get_next_city(vector<double> &probabilities);
    void calculate_ant_routes(Ant *ant, vector<vector<int>> &routes, vector<vector<double>> &pheromones);
    void pa();
    void ga();
    };


#endif //PEA_0_GRAPH_H