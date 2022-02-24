#include <fstream>
#include <cmath>
#include "Algorithm.hpp"
#include "../utils/constants.hpp"
#include <ostream>
#include <iostream>
#include <climits>
#include <random>


Algorithm::Algorithm() {
    dimension = 0;
    min_cost = INT_MAX;
}

int Algorithm::loadData(ProblemType version,
                        string file_name) {
    fstream file;
    string line;

    switch (version) {
        case TSP:
            file_name = string(TSP_PATH) + file_name;
            break;
        case SMALL:
            file_name = string(SMALL_PATH) + file_name;
            break;
        case ATSP:
            file_name = string(ATSP_PATH) + file_name;
            break;
        default:
            return EXEC_ERROR;
    }
    file.open(file_name, ios::in);
    if (!file.good()) {
        cout << "\nCouldn't open a file!\n";
        file.close();
        file.clear();
        return EXEC_ERROR;
    } else {
        getline(file, line);
        if (line.find(".tsp") == string::npos && version == TSP) {
            cout << "TSP file not found!" << endl;
            return EXEC_ERROR;
        } else if (line.find(".atsp") == string::npos && version == ATSP) {
            cout << "ATSP file not found!" << endl;
            return EXEC_ERROR;
        } else if (line.find("data") == string::npos && version == SMALL) {
            cout << "SMALL file not found!" << endl;
            return EXEC_ERROR;
        }
        getline(file, line);
        dimension = stoi(line);

        matrix = vector<vector<int>>(dimension);
        for (int i = 0; i < dimension; i++) {
            matrix[i] = vector<int>(dimension);
        }

        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                file >> matrix[i][j];
            }
        }

        file.close();
        file.clear();
        return 0;
    }
}

int Algorithm::printGraph() {
    cout << "Loaded graph: " << endl;
    if (matrix.empty()) {
        cout << "Nothing to print!\n";
        return EXEC_ERROR;
    }
    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            cout.width(PRINTED_FIELD_SIZE);
            cout << matrix[i][j];
        }
        cout << endl;
    }
    return 0;
}

int Algorithm::calculateCost(const vector<int> &vec) {
    int cost = 0;
    if (matrix.empty() || vec.empty()) {
        cout << "Nothing to count!\n";
        return EXEC_ERROR;
    }
    for (int i = 0; i < dimension - 1; i++) {
        cost += matrix[vec[i]][vec[i + 1]];
    }
    cost += matrix[vec[dimension - 1]][vec[0]];

    return cost;
}

int Algorithm::getDimension() const {
    return dimension;
}

static inline int swap_linear_update(int &current_val, vector<int> &vec,
                                     int &i, int &j, vector<vector<int>> &matrix) {
    return current_val - matrix[vec[i - 1]][vec[i]]
           - matrix[vec[i]][vec[i + 1]]
           - matrix[vec[j - 1]][vec[j]]
           - matrix[vec[j]][vec[j + 1]]
           + matrix[vec[i - 1]][vec[j]]
           + matrix[vec[j]][vec[i + 1]]
           + matrix[vec[j - 1]][vec[i]]
           + matrix[vec[i]][vec[j + 1]];
}

static inline int invert_linear_update(int &current_val, vector<int> &vec,
                                       int &i, int &j, vector<vector<int>> &matrix) {
    return current_val - matrix[vec[i - 1]][vec[i]]
           - matrix[vec[j]][vec[j + 1]]
           + matrix[vec[i - 1]][vec[j]]
           + matrix[vec[i]][vec[j + 1]];
}

static inline int insert_linear_update(int &current_val, vector<int> &vec,
                                       int &i, int &j, vector<vector<int>> &matrix) {
    return current_val - matrix[vec[i - 1]][vec[i]]
           - matrix[vec[j - 1]][vec[j]]
           - matrix[vec[j]][vec[j + 1]]
           + matrix[vec[j]][vec[i]]
           + matrix[vec[i - 1]][vec[j]]
           + matrix[vec[j - 1]][vec[j + 1]];
}

int Algorithm::greedy(vector<int> &arr, int starting_point) {
    bool visited_cities[dimension];
    for (int i = 0; i < dimension; i++) visited_cities[i] = false;

    int cost = 0;
    int aux_cost;
    int best_opt;
    visited_cities[starting_point] = true;
    int current_pos = starting_point;
    arr[0] = starting_point;

    for (int i = 0; i < dimension; i++) {
        if (i == dimension - 1) {
            cost += matrix[current_pos][starting_point];
            return cost;
        }

        aux_cost = INT_MAX;

        for (int j = 0; j < dimension; j++) {
            if (current_pos == j || visited_cities[j]) continue;
            if (matrix[current_pos][j] < aux_cost) {
                aux_cost = matrix[current_pos][j];
                best_opt = j;
            }
        }
        visited_cities[best_opt] = true;
        arr[i + 1] = best_opt;
        current_pos = best_opt;
        cost += aux_cost;
    }
    return cost;
}

int Algorithm::genRandomPathCost(vector<int> &vec) {
    shuffle(begin(vec), end(vec), std::mt19937(std::random_device()()));
    return calculateCost(vec);
}

void Algorithm::chooseCommonParams(vector<int> &vec, NeighborhoodType &nt, BeginningSolution &result, move_foo &move,
                                   move_linear_update &linear_update) {
    char option;

    cout << "\n";
    cout << "Choose beginning solution (random by default): \n";
    cout << " 1 - Natural.\n";
    cout << " 2 - Random.\n";
    cout << " 3 - Greedy.\n";
    cout << "\n";
    cin >> option;

    switch (option) {
        case '1':
            result = NATURAL;
            break;
        case '2':
            result = RANDOM;
            break;
        case '3':
            result = GREEDY; //greedy(arr, (m_rand() % (dimension + 1)));
            break;
        default:
            result = RANDOM;
            break;
    }

    cout << "\n";
    cout << "Choose neighbourhood type (swap by default): \n";
    cout << " 1 - Swap.\n";
    cout << " 2 - Invert.\n";
    cout << " 3 - Insert.\n";
    cout << "\n";
    cin >> option;

    switch (option) {
        case '1':
            move = swap;
            linear_update = swap_linear_update;
            nt = SWAP;
            break;
        case '2':
            move = invert;
            linear_update = invert_linear_update;
            nt = INVERT;
            break;
        case '3':
            move = insert;
            linear_update = insert_linear_update;
            nt = INSERT;
            break;
        default:
            move = swap;
            break;
    }
}

int Algorithm::getNextCity(vector<double> &probabilities) {
    double val = (double) m_rand(0, RAND_MAX) / ((double) (RAND_MAX + 1.0));

    int i = 0;
    double sum = probabilities[i];
    while (sum < val) { // choose next city where the sum of probabilities meets the random value
        i++;
        sum += probabilities[i];
    }

    return i;
}

void Algorithm::printBestPath() {
    cout << "Best path: ";
    for (int i = 0; i < dimension; i++) {
        if (i < dimension - 1) {
            cout << best_path[i] << ", ";
        } else {
            cout << best_path[i] << "\n";
        }
    }
}

void Algorithm::printDuration(chrono::steady_clock::time_point begin, chrono::steady_clock::time_point end) {
    cout << "Duration time: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
         << "[ms]" << "\n\n";
}


void Algorithm::printSummary(chrono::steady_clock::time_point begin, chrono::steady_clock::time_point end) {
    printBestPath();
    printDuration(begin, end);
}
