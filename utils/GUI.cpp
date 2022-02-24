//
// Created by Mariusz Wiśniewski on 2/23/2022.
//

#include "GUI.hpp"
#include "../algorithms/Algorithm.hpp"
#include "../algorithms/brute_force/SwapBruteForce.hpp"
#include "../algorithms/brute_force/TreeSearchBruteForce.hpp"
#include "../algorithms/dynamic_programming/DynamicPrograming.hpp"
#include "../algorithms/branch_and_bound/BranchAndBound.hpp"
#include "../algorithms/tabu_search/TabuSearch.hpp"
#include "../algorithms/simulated_annealing/SimulatedAnnealing.hpp"
#include "../algorithms/ant_colony_optimization/AntColonyOptimization.hpp"
#include "../algorithms/genetic_algorithm/GeneticAlgorithm.hpp"

#include <vector>
#include <iostream>

int GUI::run() {
    char option;
    string temp_opt;
    int rc = 0;
    vector<int> vec;
    string file_name;
    static Algorithm algorithm;

    temp_opt = displayMenuStart();

    if (temp_opt.size() != 1) {
        option = EXEC_ERROR;
    } else {
        option = temp_opt[0];
    }

    switch (option) {
        case '1':
            cout << "\n";
            cout << "Choose problem's version: \n";
            cout << " 1 - TSP.\n";
            cout << " 2 - ATSP.\n";
            cout << " 3 - SMALL.\n";
            cout << " 4 - Go back to main menu.\n";
            cout << "\n";
            do {
                cout << " Introduce your choice and press Enter: ";
                cin >> temp_opt;
                fflush(stdin);
                if (temp_opt.size() != 1) {
                    option = EXEC_ERROR;
                } else {
                    option = temp_opt[0];
                }
                if (option == '1' || option == '2' || option == '3') {
                    cout << "Please enter a valid file name for chosen problem's version: ";
                    cin >> file_name;
                    cout << "\n";
                }

                ProblemType type;
                switch (option) {
                    case '1':
                        type = TSP;
                        break;
                    case '2':
                        type = ATSP;
                        break;
                    case '3':
                        type = SMALL;
                        break;
                    case '4':
                        return 0;
                    default:
                        cout << "Wrong option, please try again\n\n";
                        rc = REPEAT;
                }

                if (rc != REPEAT) {
                    rc = algorithm.loadData(type, file_name);
                    if (rc == EXEC_ERROR) {
                        cout << "Loading error, please verify your input data and try again.\n";
                        return 0;
                    } else cout << "Data loaded correctly!\n";
                }
            } while (rc == REPEAT);

            break;
        case '2':
            rc = algorithm.printGraph();

            if (rc == EXEC_ERROR) {
                cout << "Printing error, please verify your input data and try again.\n";
                return 0;
            }
            break;
        case '3':
            vec = vector<int>(algorithm.getDimension());
            for (int i = 0; i < algorithm.getDimension(); i++) {
                vec[i] = i;
            }
            rc = algorithm.calculateCost(vec);

            if (rc == EXEC_ERROR) {
                cout
                        << "Error during displaying the cost, please verify if your input data is correct and try again.\n";
                return 0;
            }
            cout << "Cost: " << rc << "\n";
            break;
        case '4': {
            SwapBruteForce sbf(algorithm);
            sbf.run();
            break;
        }
        case '5' : {
            TreeSearchBruteForce tsbf(algorithm);
            tsbf.run();
            break;
        }
        case '6': {
            DynamicProgramming dp(algorithm);
            dp.run();
            break;
        }
        case '7': {
            BranchAndBound bab(algorithm);
            bab.run();
            break;
        }
        case '8': {
            TabuSearch ts(algorithm);
            ts.run();
            break;
        }
        case '9': {
            SimulatedAnnealing sa(algorithm);
            sa.run();
            break;
        }
        case '0':
            return EXIT;
        case 'a': {
            char option;

            cout << "\n";
            cout << "Choose pheromones update type (density by default): \n";
            cout << " 1 - DAS.\n";
            cout << " 2 - CAS.\n";
            cout << " 3 - QAS.\n";
            cout << "\n";
            cin >> option;

            AntColonyOptimizationType type;
            switch (option) {
                case '1': {
                    type = AntColonyOptimizationType::DENSITY;
                    break;
                }
                case '2': {
                    type = AntColonyOptimizationType::CYCLE;
                    break;
                }
                case '3': {
                    type = AntColonyOptimizationType::QUANTITY;
                    break;
                }
                default: {
                    type = AntColonyOptimizationType::DENSITY;
                    break;
                }
            }

            AntColonyOptimization aco(algorithm, type);
            aco.run();
            break;
        }
        case 'b': {
            GeneticAlgorithm ga(algorithm);
            ga.run();
            break;
        }
        default:
            cout << "Wrong option, try again!\n";
            break;
    }
    return 0;
}

string GUI::displayMenuStart() {
    string opt;

    cout << "\n";
    cout << "****************************************\n";
    cout << " 1 - Load a graph.\n";
    cout << " 2 - Display loaded graph.\n";
    cout << " 3 - Display the cost.\n";
    cout << " 4 - Run the brute force based on swaps.\n";
    cout << " 5 - Run the brute force based on tree.\n";
    cout << " 6 - Run the dynamic programming.\n";
    cout << " 7 - Run BnB.\n";
    cout << " 8 - Run Tabu Search.\n";
    cout << " 9 - Run Simulated Annealing.\n";
    cout << " a - Run Ant Colony Algorithm.\n";
    cout << " b - Run Genetic Algorithm.\n";
    cout << " 0 - Exit.\n";
    cout << " Introduce your choice and press Enter: ";
    cin >> opt;
    fflush(stdin);

    return opt;
}
