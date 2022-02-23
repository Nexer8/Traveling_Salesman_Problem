cmake_minimum_required(VERSION 3.15)
project(TSP)

set(CMAKE_CXX_STANDARD 17)

add_executable(TSP
        main.cpp
        Algorithm.cpp
        Algorithm.hpp
        utils/constants.hpp
        algorithms/ant_colony_optimization/Ant.cpp
        algorithms/ant_colony_optimization/Ant.hpp
        algorithms/simulated_annealing/SimulatedAnnealing.cpp
        algorithms/simulated_annealing/SimulatedAnnealing.hpp
        algorithms/tabu_search/TabuSearch.cpp
        algorithms/tabu_search/TabuSearch.hpp
        algorithms/branch_and_bound/BranchAndBound.cpp
        algorithms/branch_and_bound/BranchAndBound.hpp
        algorithms/ant_colony_optimization/AntColonyOptimization.cpp
        algorithms/ant_colony_optimization/AntColonyOptimization.hpp
        algorithms/genetic_algorithm/GeneticAlgorithm.cpp
        algorithms/genetic_algorithm/GeneticAlgorithm.hpp
        algorithms/brute_force/SwapBruteForce.cpp
        algorithms/brute_force/SwapBruteForce.hpp
        algorithms/brute_force/TreeSearchBruteForce.cpp
        algorithms/brute_force/TreeSearchBruteForce.hpp
        algorithms/dynamic_programming/DynamicPrograming.cpp
        algorithms/dynamic_programming/DynamicPrograming.hpp
        utils/GUI.cpp
        utils/GUI.hpp
)