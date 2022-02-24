#include <iostream>
#include <string>
#include <vector>
#include <functional>
#include <random>
#include <chrono>
#include "../utils/constants.hpp"

#ifndef PEA_0_GRAPH_H
#define PEA_0_GRAPH_H

#define EXEC_ERROR -1

using namespace std;


typedef void (*crossover_type)(vector<int> &, vector<int> &);

class Algorithm {
protected:
    int dimension;
    vector<vector<int>> matrix;
    int min_cost;
    vector<int> best_path;
    vector<vector<int>> population;
    vector<int> fitness;

    static inline void swap(vector<int> &vec, int a, int b) {
        int temp = vec[a];
        vec[a] = vec[b];
        vec[b] = temp;
    }

    static inline void insert(vector<int> &vec, int dst, int src) {
        int element_to_insert = vec.at(src);
        if (dst > src) {
            for (unsigned int i = src; i < dst; i++) {
                vec.at(i) = vec.at(i + 1);
            }
        } else {
            for (unsigned int i = src; i > dst; i--) {
                vec.at(i) = vec.at(i - 1);
            }
        }
        vec.at(dst) = element_to_insert;
    }

    static inline void invert(vector<int> &vec, int low, int high) {
        if (low > high) {
            swap(vec, low, high);
        }
        for (int i = low, j = high; i < j; i++, j--) {
            swap(vec, i, j);
        }
    }

    static void chooseCommonParams(vector<int> &vec, NeighborhoodType &nt, BeginningSolution &result, move_foo &move,
                                   move_linear_update &linear_update);

    int genRandomPathCost(vector<int> &vec);

    int greedy(vector<int> &arr, int starting_point);

    int getNextCity(vector<double> &probabilities);

public:
    int calculateCost(const vector<int> &vec);

    int getDimension() const;

    int loadData(ProblemType version, string file_name);

    int printGraph();

    Algorithm();

    virtual void run() {};

    void printBestPath();

    static void printDuration(chrono::steady_clock::time_point begin, chrono::steady_clock::time_point end);

    void printSummary(chrono::steady_clock::time_point begin, chrono::steady_clock::time_point end);

    static inline std::mt19937 &generator() {
        static thread_local std::mt19937 gen(std::random_device{}());
        return gen;
    }

    template<typename T, std::enable_if_t<std::is_integral_v<T>> * = nullptr>
    static T m_rand(T min, T max) {
        std::uniform_int_distribution<T> dist(min, max);
        return dist(generator());
    }

    template<typename T, std::enable_if_t<std::is_floating_point_v<T>> * = nullptr>
    static T m_rand(T min, T max) {
        std::uniform_real_distribution<T> dist(min, max);
        return dist(generator());
    }
};

#endif //PEA_0_GRAPH_H
