#include <chrono>
#include "Graph.h"

#define EXIT -2
#define REPEAT 1
using namespace std;

int print_menu(Graph *cgraph) {
    char option;
    string temp_opt;
    int rc;
    int *arr;
    string file_name;
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
    cin >> temp_opt;
    fflush(stdin);
    if (temp_opt.size() != 1) {
        option = EXEC_ERROR;
    }
    else {
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
                }
                else {
                    option = temp_opt[0];
                }
                if (option == '1' || option == '2' || option == '3') {
                    cout << "Please enter a valid file name for chosen problem's version: ";
                    cin >> file_name;
                    cout << "\n";
                }
                switch (option) {
                    case '1':
                        rc = cgraph->load_data(TSP, file_name);
                        if (rc == EXEC_ERROR) {
                            cout << "Loading error, please verify your input data and try again.\n";
                            return 0;
                        }
                        else cout << "Data loaded correctly!\n";
                        break;
                    case '2':
                        rc = cgraph->load_data(ATSP, file_name);
                        if (rc == EXEC_ERROR) {
                            cout << "Loading error, please verify your input data and try again.\n";
                            return 0;
                        }
                        else cout << "Data loaded correctly!\n";
                        break;
                    case '3':
                        rc = cgraph->load_data(SMALL, file_name);
                        if (rc == EXEC_ERROR) {
                            cout << "Loading error, please verify your input data and try again.\n";
                            return 0;
                        }
                        else cout << "Data loaded correctly!\n";
                        break;
                    case '4':
                        return 0;
                    default:
                        cout << "Wrong option, please try again\n\n";
                        rc = REPEAT;
                }
            } while (rc == REPEAT);

            break;
        case '2':
            rc = cgraph->print_graph();

            if (rc == EXEC_ERROR) {
                cout << "Printing error, please verify your input data and try again.\n";
                return 0;
            }
            break;
        case '3':
            /* TODO: This part needs to be moved to another function that will generate permutations */
            arr = (int *)malloc(sizeof(int) * cgraph->get_dimension());
            for (int i = 0; i < cgraph->get_dimension(); i++) {
                arr[i] = i;
            }
            rc = cgraph->calculate_cost(arr);
            free(arr);

            if (rc == EXEC_ERROR) {
                cout << "Error during displaying the cost, please verify if your input data is correct and try again.\n";
                return 0;
            }
            cout << "Cost: " << rc << "\n";
            break;
        case '4': {
            if (cgraph->get_dimension() != 0) {
                arr = (int *) malloc(cgraph->get_dimension() * sizeof(int));
                cgraph->best_path = (int *) malloc(cgraph->get_dimension() * sizeof(int));
                for (int i = 0; i < cgraph->get_dimension(); i++) {
                    arr[i] = i;
                }
                cgraph->min_cost = INT_MAX;
                chrono::steady_clock::time_point begin_bf = chrono::steady_clock::now();
                cgraph->brute_force_swap(0, arr);
                chrono::steady_clock::time_point end_bf = chrono::steady_clock::now();

                cout << "Minimal cost: " << cgraph->min_cost << "\n";
                cout << "Best path: ";
                for (int i = 0; i < cgraph->get_dimension(); i++) {
                    if (i != cgraph->get_dimension() - 1) {
                        cout << cgraph->best_path[i] << ", ";
                    } else {
                        cout << cgraph->best_path[i];
                    }
                }
                cout << "\nDuration time: "
                     << chrono::duration_cast<std::chrono::milliseconds>(end_bf - begin_bf).count() << "[ms]" << "\n";

                cout << "\n";
                free(arr);
                free(cgraph->best_path);
            } else {
                cout << "Error while trying to read the matrix.\n";
            }
            break;
        }
        case '5' : {
            cgraph->best_path = (int *) malloc(cgraph->get_dimension() * sizeof(int));
            chrono::steady_clock::time_point begin = chrono::steady_clock::now();
            cgraph->brute_force_tree();
            chrono::steady_clock::time_point end = chrono::steady_clock::now();

            cout << "Best path: ";
            for (int i = 0; i < cgraph->get_dimension(); i++) {
                if (i != cgraph->get_dimension() - 1) {
                    cout << cgraph->best_path[i] << ", ";
                } else {
                    cout << cgraph->best_path[i];
                }
            }

            cout << "\nDuration time: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
                 << "[ms]" << "\n";
            cout << "\n";
            free(cgraph->best_path);
            break;
        }
        case '6': {
            cgraph->best_path = (int *) malloc(cgraph->get_dimension() * sizeof(int));
            cgraph->min_cost = INT_MAX;

            chrono::steady_clock::time_point begin_dp = chrono::steady_clock::now();
            cgraph->dynamic_programming();
            chrono::steady_clock::time_point end_dp = chrono::steady_clock::now();

            cout << "\nDuration time: " << chrono::duration_cast<std::chrono::milliseconds>(end_dp - begin_dp).count()
                 << "[ms]" << "\n";
            cout << "\n";
            free(cgraph->best_path);
            break;
        }
        case '7':
            cgraph->branch_n_bound();
            break;
        case '8': {
            cgraph->tabu_search();
            break;
        }
        case '9': {
            cgraph->simulated_annealing();
            break;
        }
        case '0':
            return EXIT;
        case 'a':
            cgraph->pa();
            break;
        case 'b':
            cgraph->ga();
            break;
        default:
            cout << "Wrong option, try again!\n";
            break;
    }
    return 0;
}

int main() {
    Graph cgraph;
    int rc = 0;
    cout << "****************************************\n";
    cout << "PEA - etap 0, Mariusz Wisniewski, 241393\n";

    while (rc != EXIT) {
        rc = print_menu(&cgraph);
    }

    return 0;
}