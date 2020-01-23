#include <fstream>
#include <cmath>
#include "Graph.h"
#include "main_program.h"
#include <chrono>
#include <c++/4.8.3/algorithm>
#include <ostream>
#include <iostream>
#include <c++/4.8.3/map>

#define STARTING_TEMPERATURE 1e9
#define STOPPING_TEMPERATURE 0.01
#define NUMBER_OF_ITERATIONS 500
#define LINEAR_COOLING_CONST 0.95
#define GEOMETRIC_COOLING_CONST 5
#define REPEAT_VAL 3

#define POPULATION_SIZE 5000
#define NUMBER_OF_GENERATIONS 500
#define NUMBER_OF_TOURNAMENTS 5
#define CROSS_RATE 0.2
#define MUTATION_RATE 0.02

#define STEPS_IN_GEN dimension
#define CADENCE dimension

Graph::Graph() {
    dimension = 0;
    best_path = nullptr;
    matrix = nullptr;
    min_cost = INT_MAX;
}

Graph::~Graph() {
    if (matrix != nullptr) {
        for (int i = 0; i < dimension; i++) {
            delete[] matrix[i];
        }
        delete[] matrix;
    }
}

int Graph::load_data(problem_type version,
              string file_name) {
    fstream file;
    string line;

    switch (version) {
        case TSP:
            file_name = current_dir + tsp_path + file_name;
            break;
        case SMALL:
            file_name = current_dir + small_path + file_name;
            break;
        case ATSP:
            file_name = current_dir + atsp_path + file_name;
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
        }
        else if (line.find(".atsp") == string::npos && version == ATSP) {
            cout << "ATSP file not found!" << endl;
            return EXEC_ERROR;
        }
        else if (line.find("data") == string::npos && version == SMALL) {
            cout << "SMALL file not found!" << endl;
            return EXEC_ERROR;
        }
        getline(file, line);
        dimension = stoi(line);

        matrix = new int *[dimension];
        for (int i = 0; i < dimension; i++) {
            matrix[i] = new int[dimension];
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

int Graph::print_graph() {
    cout << "Loaded graph: " << endl;
    if (matrix == nullptr) {
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

int Graph::calculate_cost(const int *arr) {
    int cost = 0;
    if (matrix == nullptr || arr == nullptr) {
        cout << "Nothing to count!\n";
        return EXEC_ERROR;
    }
    for (int i = 0; i < dimension - 1; i++) {
        cost += matrix[arr[i]][arr[i + 1]];
    }
    cost += matrix[arr[dimension - 1]][arr[0]];

    return cost;
}

int Graph::get_dimension() {
    return dimension;
}

static inline void swap(int *arr, int a, int b) {
    int temp = arr[a];
    arr[a] = arr[b];
    arr[b] = temp;
}

static inline void insert(int *arr, int dst, int src) {
    int dim  = abs(src - dst);
    int temp[dim];
    int idx = 0;

    if (src > dst) {
        for (int i = dst; i <= src; i++) {
            temp[idx] = arr[i];
            idx++;
        }
        arr[dst] = arr[src];
        idx = 0;
        for (int i = dst + 1; i <= src; i++) {
            arr[i] = temp[idx];
            idx++;
        }
    }
    else if (src < dst) {
        idx = 0;
        for (int i = src; i <= dst; i++) {
            temp[idx] = arr[i];
            idx++;
        }
        arr[dst] = arr[src];
        idx = 0;
        for (int i = src; i < dst; i++) {
            arr[i] = temp[idx + 1];
            idx++;
        }
    }
    else { return; }
}

static inline void invert(int *arr, int low, int high) {
    if (low > high) {
        swap(arr, low, high);
    }
    for (int i = low, j = high; i < j; i++, j--) {
        swap(arr, i, j);
    }
}

static inline int swap_linear_update(int* current_val, int* arr, int* i, int* j, int **matrix) {
    return *current_val - matrix[arr[*i - 1]][arr[*i]]
           - matrix[arr[*i]][arr[*i + 1]]
           - matrix[arr[*j - 1]][arr[*j]]
           - matrix[arr[*j]][arr[*j + 1]]
           + matrix[arr[*i - 1]][arr[*j]]
           + matrix[arr[*j]][arr[*i + 1]]
           + matrix[arr[*j - 1]][arr[*i]]
           + matrix[arr[*i]][arr[*j + 1]];
}

static inline int invert_linear_update(int* current_val, int* arr, int* i, int* j, int **matrix) {
    return *current_val - matrix[arr[*i - 1]][arr[*i]]
           - matrix[arr[*j]][arr[*j + 1]]
           + matrix[arr[*i - 1]][arr[*j]]
           + matrix[arr[*i]][arr[*j + 1]];
}

static inline int insert_linear_update(int* current_val, int* arr, int* i, int* j, int **matrix) {
    return *current_val - matrix[arr[*i - 1]][arr[*i]]
                        - matrix[arr[*j - 1]][arr[*j]]
                        - matrix[arr[*j]][arr[*j + 1]]
                        + matrix[arr[*j]][arr[*i]]
                        + matrix[arr[*i - 1]][arr[*j]]
                        + matrix[arr[*j - 1]][arr[*j + 1]];
}

void Graph::brute_force_swap(int current_level, int *arr) {
    if (matrix == nullptr || current_level < 0 || current_level > dimension) return;
    if (current_level == dimension) {
        int cost = calculate_cost(arr);
        if (cost < min_cost) {
            min_cost = cost;
            for (int i = 0; i < dimension; i++) {
                best_path[i] = arr[i];
            }
        }
    }

    for (int i = current_level; i < dimension; i++) {
        // We swap the elements to generate a arr
        swap(arr, current_level, i);
        brute_force_swap(current_level + 1, arr);
        // Then we swap them back - backtracking
        swap(arr, current_level, i);
    }
}

void Graph::tree_search(vector<int> vertex_vector) {
    if (vertex_vector.size() == dimension) { // when the path is complete we check its cost
        int arr[dimension];
        int k = 0;
        for (int const &i: vertex_vector) {
            arr[k++] = i;
        }
        int temp_result = calculate_cost(arr);
        if (temp_result < min_cost) {
            min_cost = temp_result;
            for (int i = 0; i < dimension; i++) {
                best_path[i] = arr[i];
            }
            return;
        }
    }

    for (int i = 0; i < dimension; i++) { // checking which vertex we can add
        bool exists = false;
        for (int j: vertex_vector) {
            if (i == j) {
                exists = true;
                break;
            }
        }
        if (!exists) { // if we can add a vertex we do this and continue the recursion
            vertex_vector.push_back(i);
            tree_search(vertex_vector);
            vertex_vector.pop_back(); // then we take the vertex out and check a new path
        }
    }
}

void Graph::brute_force_tree() {
    vector<int> vertex_vector;
    min_cost = INT_MAX;
    vertex_vector.push_back(0);
    tree_search(vertex_vector);
    cout << "Minimal cost: " << min_cost << "\n";
}

void Graph::dynamic_programming() {
    int **dp_matrix = nullptr; // [path][last vertex]
    dp_matrix = new int*[1u << (unsigned)dimension]; // create a room for all possible paths

    for (unsigned int i = 0; i < (1u << (unsigned)dimension); i++) {
        dp_matrix[i] = new int[dimension]; // number of columns is equal to matrix dimension
    }

    for (unsigned int i = 0; i < (1u << (unsigned)dimension); i++) {
        for (int j = 0; j < dimension; j++) {
            dp_matrix[i][j] = 1u << 30u; // set all the costs to "infinity"
        }
    }

    for (unsigned int i = 0; i < dimension; i++) {
        if (i == 0) dp_matrix[1u << i | 1u][i] = 0; // Don't want a -1 value to be put in here.
        else {
            dp_matrix[1u << i | 1u][i] = matrix[0][i]; // a cost of 0<->(n-1) path
        }
    }

    for (unsigned int bit_mask = 0; bit_mask < 1u << (unsigned)dimension; bit_mask++) {
        for (unsigned int vertex = 0; vertex < (unsigned)dimension; vertex++) {
            if (!(bit_mask & (1u << vertex))) continue; // filter just the paths that don't include the vertex and has to end in it
            for (unsigned int j = 0; j < (unsigned)dimension; j++) {
                if (!(bit_mask & (1u << j))) {// vertex j which is a proposed next extension cannot already exist in the path
                    if (dp_matrix[bit_mask][vertex] + matrix[vertex][j] < dp_matrix[bit_mask | (1u << j)][j]) { // check if it's better to add vertex j to the current path that finishes in vertex or the actual path till j is more optimal
                        dp_matrix[bit_mask | (1u << j)][j] = dp_matrix[bit_mask][vertex] + matrix[vertex][j];
                    }
                }
            }
        }
    }

    int cost;

    for (int i = 0; i < dimension; i++) { // looking for the minimal cost cycle in the last row
        cost = dp_matrix[(1u << (unsigned)dimension) - 1][i] + matrix[i][0];
        if (cost < min_cost) {
            min_cost = cost;
        }
    }

    cout << "Minimal cost: " << min_cost << "\n";
}

int Graph::min_distance_from(int vertex) {
    int min = INT_MAX;
    for (int i = 0; i < dimension; i++)
    {
        if (matrix[vertex][i] < min && i != vertex) min = matrix[vertex][i]; // finding the minimum edge cost having an end at the vertex "vertex"
    }
    return min;
}

int Graph::second_distance_from(int vertex) {
    int min = INT_MAX, second = INT_MAX;
    for (int i = 0; i < dimension; i++)
    {
        if (i == vertex) continue;

        if (matrix[vertex][i] < min && i != vertex) {
            min = matrix[vertex][i]; // finding the second minimum edge cost having an end at the vertex "vertex"
            second = min;
        }
        else if (matrix[vertex][i] <= second && matrix[vertex][i] != min) second = matrix[vertex][i];
    }
    return second;
}

void Graph::branch_n_bound() {
    int *tab = new int[dimension];
    for (int i = 0; i < dimension; i++) tab[i] = i;

    int lower_bound = 0;
    // setting the lower bound so that is the sum of minimal distances for each vertex
    // We divide it by two because if not then we go through each path twice
    for (int i = 0; i < dimension; i++) lower_bound += ceil((min_distance_from(i) + second_distance_from(i)) / 2.0);

    int level = 1;
    int distance = 0;
    min_cost = INT_MAX;
    best_path = new int[dimension];
    chrono::steady_clock::time_point begin = chrono::steady_clock::now();

    branch_n_bound(tab, level, lower_bound, distance, min_cost);

    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    cout << "Minimal cost: " << min_cost << "\n";
    cout << "Best path: ";
    for (int i = 0; i < dimension; i++) {
        if (i < dimension - 1) {
            cout << best_path[i] << ", ";
        }
        else {
            cout << best_path[i] << "\n";
        }
    }

    cout << "Duration time: " << chrono::duration_cast<chrono::milliseconds>(end - begin).count() << "[ms]" << "\n";

    cout << "\n";

    delete[] tab;
    delete[] best_path;
}

void Graph::branch_n_bound(int *tab, int level, int lower_bound, int distance, int& min_distance) {
    if (level == dimension) {
        if (matrix[tab[level - 1]][tab[0]] < 0) return;

        // checking if the distance is smaller at the end
        distance += matrix[tab[level - 1]][tab[0]];
        if (distance < min_distance) {
            for (int i = 0; i < dimension; i++) best_path[i] = tab[i];
            min_distance = distance; // upper bound
        }
    }
    else {
        for (int i = level; i < dimension; i++) {
            if (matrix[tab[level - 1]][tab[i]] < 0) continue;

            // adjusting the lower_bound so that it fits the level we're currently at
            int temp_bound = lower_bound;
            if (level == 1) lower_bound -= ceil((min_distance_from(tab[0]) + min_distance_from(tab[i])) / 2.0);
            else if (level == dimension - 1) lower_bound -= ceil((second_distance_from(tab[level - 1]) + second_distance_from(tab[0])) / 2.0);
            else lower_bound -= ceil((second_distance_from(tab[level - 1]) + min_distance_from(tab[i])) / 2.0);

            // distance from 0 to i
            distance += matrix[tab[level - 1]][tab[i]];

            // verify if further checking makes sense
            // lower_bound + actual distance should be then smaller than the minimal distance found, if not, we cut the node
            if (lower_bound + distance < min_distance) {
                swap(tab, level, i);

                branch_n_bound(tab, level + 1, lower_bound, distance, min_distance);

                swap(tab, level, i);
            }

            distance -= matrix[tab[level - 1]][tab[i]]; // we go back in a tree
            lower_bound = temp_bound;
        }
    }
}

int Graph::greedy(int *arr, int starting_point) {
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

int Graph::gen_random_path_cost(int *arr) {
    random_shuffle(&arr[0], &arr[dimension]);
    return calculate_cost(arr);
}

void Graph::ts_choose_params(int *arr, neighbourhood_type *nt, beginning_solution *result, move_foo *move, move_linear_update *linear_update) {
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
            *result = NATURAL;
            break;
        case '2':
            *result = RANDOM;
            break;
        case '3':
            *result = GREEDY; //greedy(arr, (rand() % (dimension + 1)));
            break;
        default:
            *result = RANDOM;
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
            *move = swap;
            *linear_update = swap_linear_update;
            *nt = SWAP;
            break;
        case '2':
            *move = invert;
            *linear_update = invert_linear_update;
            *nt = INVERT;
            break;
        case '3':
            *move = insert;
            *linear_update = insert_linear_update;
            *nt = INSERT;
            break;
        default:
            *move = swap;
            break;
    }
}

void Graph::tabu_search() {
    move_foo move;
    move_linear_update linear_update;
    neighbourhood_type ngbh_type;

    best_path = new int[dimension];
    beginning_solution beg_solution;
    int result;
    int arr[dimension];
    for (int i = 0; i < dimension; i++) {
        arr[i] = i;
    }

    // Setting a type of the neighbourhood and the beginning solution
    ts_choose_params(arr, &ngbh_type, &beg_solution, &move, &linear_update);
    chrono::steady_clock::time_point begin = chrono::steady_clock::now();

    switch (beg_solution) {
        case NATURAL:
            result = calculate_cost(arr);
            break;
        case RANDOM:
            result = gen_random_path_cost(arr);
            break;
        case GREEDY:
            result = greedy(arr, 0);
            break;
    }
    int current_val = result;
    int temp_path[dimension];

    // Creating a tabu_list
    vector<vector<int>> tabu_list;
    tabu_list.resize(dimension);
    for (vector<vector<int> >::iterator it = tabu_list.begin(); it != tabu_list.end(); it++) {
        it->resize(dimension, 0);
    }

    int next_step[dimension];

    for (int num_of_iteration = 0; num_of_iteration < NUMBER_OF_ITERATIONS; num_of_iteration++) { // setting a number of iterations
        for (int step = 0; step < STEPS_IN_GEN * dimension; step++) {   // number of steps in a generation
            int next_step_val = INT_MAX;

            int f_tabu = 0, s_tabu = 0;

            for (int i = 0; i < dimension; i++) {
                for (int j = i + 1; j < dimension; j++) {
                    int current_val_tmp = current_val;

                    if (ngbh_type == INSERT) {
                        memcpy(temp_path, arr, dimension*sizeof(int)); // copy solution array to best_path array
                    }
                    if ((j - i < 3) || (j == dimension - 1) || (i == 0) && ngbh_type != INSERT) { // boundary cases, just update in O(n) time
                        move(arr, i, j);
                        current_val = calculate_cost(arr);
                    }
                    else if ((i == 0 || j == dimension - 1) && ngbh_type == INSERT) {
                        move(arr, i, j);
                        current_val = calculate_cost(arr);
                    }
                    else {  // O(1) update
                        current_val = linear_update(&current_val, arr, &i, &j, matrix);
                        move(arr, i, j);
                    }

                    if (current_val < result) { // if the solution in neighbourhood is better than the current global best
                        result = current_val; // aspiration criteria - don't even check the tabu_list
                        memcpy(best_path, arr, dimension*sizeof(int)); // copy solution array to best_path array
                    }

                    if (current_val < next_step_val) {  // if the solution is better than the local best (inside a given step)
                        if (tabu_list[j][i] <= step) {  // if the move is available, then update tabu_list and next_step_val
                            s_tabu = i;
                            f_tabu = j;
                            memcpy(next_step, arr, dimension*sizeof(int));
                            next_step_val = current_val;
                        }
                    }
                    if (ngbh_type != INSERT) move(arr, i, j); // go back to the previous position
                    else {
                        memcpy(arr, temp_path, dimension*sizeof(int));
                    }
                    current_val = current_val_tmp;
                }
            }
            memcpy(arr, next_step, dimension*sizeof(int));  // assigning best step
            tabu_list[f_tabu][s_tabu] += CADENCE; // updating tabu list with the opposite value
            // tabu on a given edge exists for <dimension> iterations
        }

        // generate a new set - critical event
        switch (beg_solution) {
            case NATURAL:
                for (int i = 0; i < dimension; i++) {
                    arr[i] = i;
                }
                current_val = calculate_cost(arr);
                break;
            case RANDOM:
                current_val = gen_random_path_cost(arr);
                break;
            case GREEDY:
                current_val = greedy(arr, 0);
                break;
        }
        // restoring tabu_list
        for(vector<vector<int> >::iterator it = tabu_list.begin(); it != tabu_list.end(); it++){
            it->clear();
        }
        tabu_list.clear();
        tabu_list.resize(dimension);
        for(vector<vector<int> >::iterator it = tabu_list.begin(); it != tabu_list.end(); it++){
            it->resize(dimension, 0);
        }
    }
    chrono::steady_clock::time_point end = chrono::steady_clock::now();

    cout << "Best Tabu Search result: " << result << "\n";
    cout << "Best path: ";
    for (int i = 0; i < dimension; i++) {
        if (i < dimension - 1) {
            cout << best_path[i] << ", ";
        }
        else {
            cout << best_path[i] << "\n";
        }
    }
    cout << "\nDuration time: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
         << "[ms]" << "\n";
    cout << "\n";
}

static inline double linear_cooling(double temperature, int it) {
    return temperature *= LINEAR_COOLING_CONST;
}

static inline double geometric_cooling(double temperature, int it) {
    return temperature -= GEOMETRIC_COOLING_CONST;
}

static inline double logarithmic_cooling(double temperature, int it) {
    return STARTING_TEMPERATURE/(log(it) + 1);
}

static inline double get_probability(int difference, double temperature) // checking how weak the solution is
{
    return exp(-1*difference/temperature);
}

void Graph::sa_choose_params(int *arr, neighbourhood_type *nt, beginning_solution *result, move_foo *move, cooling_method *cooling, move_linear_update *linear_update) {
    ts_choose_params(arr, nt, result, move, linear_update);
    char option;

    cout << "\n";
    cout << "Choose cooling method (linear by default): \n";
    cout << " 1 - Linear.\n";
    cout << " 2 - Geometric.\n";
    cout << " 3 - Logarithmic.\n";
    cout << "\n";
    cin >> option;

    switch (option) {
        case '1':
            *cooling = linear_cooling;
            break;
        case '2':
            *cooling = geometric_cooling;
            break;
        case '3':
            *cooling = logarithmic_cooling;
            break;
        default:
            *cooling = linear_cooling;
            break;
    }
}

void Graph::simulated_annealing() {
    move_foo move;
    cooling_method cooling;
    move_linear_update linear_update;
    neighbourhood_type ngbh_type;

    best_path = new int[dimension];
    int arr[dimension];
    for (int i = 0; i < dimension; i++) {
        arr[i] = i;
    }
    beginning_solution beg_solution;
    int result;

    sa_choose_params(arr, &ngbh_type, &beg_solution, &move, &cooling, &linear_update);

    chrono::steady_clock::time_point begin = chrono::steady_clock::now();

    switch (beg_solution) {
        case NATURAL:
            result = calculate_cost(arr);
            break;
        case RANDOM:
            result = gen_random_path_cost(arr);
            break;
        case GREEDY:
            result = greedy(arr, rand() % (dimension ));
            break;
    }

    int step_val = result;

    int next_step[dimension];
    int temp_arr[dimension];

    int stagnation = 0; // stagnation value
    int iteration = 0;  // number of iterations for one annealing and cooling
    double temperature = STARTING_TEMPERATURE;

    for (int num_of_iteration = 0; num_of_iteration < NUMBER_OF_ITERATIONS; ++num_of_iteration) {
        while (temperature > STOPPING_TEMPERATURE) { // end condition
            iteration++;
            int repeat = REPEAT_VAL*dimension;  // number of neighbourhood verifications = 3n, number of iterations for the same temperature
            int temp_step_val = step_val;

            while (repeat-- > 0) {
                int i = 0;  // first iterator
                int j = 0;  // second iterator
                while(i == j) {
                    i = rand() % dimension;
                    j = rand() % dimension;
                }
                if(i > j) swap(i, j);

                if (ngbh_type == INSERT) {
                    memcpy(temp_arr, arr, dimension*sizeof(int)); // copy solution array to best_path array
                }
                 // checking step_value, O(n) update
                if ((j - i < 3) || j == dimension - 1 || i == 0 && ngbh_type != INSERT) { // border cases verification
                    move(arr, i, j);
                    step_val = calculate_cost(arr);
                }
                else if ((i == 0 || j == dimension - 1) && ngbh_type == INSERT) {
                    move(arr, i, j);
                    step_val = calculate_cost(arr);
                } // O(1) update
                else {
                    step_val = linear_update(&step_val, arr, &i, &j, matrix);
                    move(arr, i, j);
                }

                int difference = step_val - result; // calculate difference between the best and actual results

                if (step_val < result) {  // if step_val is better than current result, update result and reset stagnation value
                    stagnation = 0;
                    result = step_val;
                    memcpy(best_path, arr, dimension*sizeof(int));
                }
                else {   // if the value has not improved, increase the stagnation value
                    stagnation++;
                    if (stagnation > dimension) temperature = (double)(STARTING_TEMPERATURE)/2;   // if stagnation is fairly big, increase the temperature
                }

                // decision making whether to take the next_step
                if (difference < 0 || (difference > 0 &&  get_probability(difference, temperature) > ( (double) rand() / (RAND_MAX)) ) ){
                    memcpy(next_step, arr, dimension*sizeof(int));
                    break;
                }
                else {
                    if (ngbh_type != INSERT) move(arr, i, j); // go back to the previous position
                    else {
                        memcpy(arr, temp_arr, dimension*sizeof(int));
                    }
                    step_val = temp_step_val;
                }
            }
            temperature = cooling(temperature, iteration);  // cooling
            memcpy(arr, next_step, dimension*sizeof(int));  // assigning best step
            if (iteration > dimension*dimension) break; // if we are for too long in the same spot
        }

        // generating new set of starting conditions
        switch (beg_solution) {
            case NATURAL:
                for (int i = 0; i < dimension; i++) {
                    arr[i] = i;
                }
                step_val = calculate_cost(arr);
                break;
            case RANDOM:
                step_val = gen_random_path_cost(arr);
                break;
            case GREEDY:
                step_val = greedy(arr, rand() % (dimension ));
                break;
        }
        iteration = 0;
        stagnation = 0;
        temperature = STARTING_TEMPERATURE;
    }
    chrono::steady_clock::time_point end = chrono::steady_clock::now();

    cout << "Best Simulated Annealing result: " << result << "\n";
    cout << "Best path: ";
    for (int i = 0; i < dimension; i++) {
        if (i < dimension - 1) {
            cout << best_path[i] << ", ";
        }
        else {
            cout << best_path[i] << "\n";
        }
    }
    cout << "\nDuration time: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
         << "[ms]" << "\n";
    cout << "\n";
}

void Graph::generate_population() {
    population.reserve(POPULATION_SIZE);
    fitness.reserve(POPULATION_SIZE);

    vector<int> route;
    for (int i = 0; i < dimension; i++) {
        route.push_back(i);
    }

    // draw new permutations and push them to the population vector
    for (int i = 0; i < POPULATION_SIZE; ++i) {
        fitness.push_back(gen_random_path_cost(&route[0]));
        population.push_back(route);
    }
}

void Graph::select_mating_pool_tournament() { // tournament mating pool selection
    vector<vector<int>> mating_pool;
    mating_pool.reserve(POPULATION_SIZE);

    for (int j = 0; j < POPULATION_SIZE; ++j) {
        int best = INT_MAX;
        int best_idx;
        for (int i = 0; i < NUMBER_OF_TOURNAMENTS; ++i) {
            int rand_idx = rand() % POPULATION_SIZE; // draw a random index
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
        fitness[i] = calculate_cost(&population[i][0]);
    }
}

static void OX(vector<int> &first_parent, vector<int> &second_parent) {
    int dimension = first_parent.size();
    int k1, k2;
    vector<int> first_child(dimension, -1);
    vector<int> second_child(dimension, -1);

    do {
        k1 = rand() % (dimension -1);
        k2 = rand() % (dimension -1);  //in case of hitting last index, while loop won't break (numbers from 1 to n-1)
    } while (k1 == k2);

    if (k1 > k2) {
        swap(k1, k2);
    }

    for (int i = k1; i <= k2; i++) {
        first_child[i] = second_parent[i];
        second_child[i] = first_parent[i];
    }

    vector<int>::iterator child_iterator = first_child.begin() + k2 + 1;
    vector<int>::iterator parent_iterator = first_parent.begin() + k2 + 1;

    while (child_iterator != first_child.begin() + k1) {
        if (first_child.end() ==
            find(first_child.begin(), first_child.end(), *parent_iterator)) { // if a child does not have parent's gene (city)
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
        } else {    // jesli miasto juz wystapilo to pomin i idz dalej
            if (parent_iterator == second_parent.end() - 1)
                parent_iterator = second_parent.begin();
            else
                parent_iterator++;
        }
    }

    first_parent = first_child;
    second_parent = second_child;
}

static void PMX(vector<int> &first_parent, vector<int> &second_parent) {
    int dimension = first_parent.size();
    int k1, k2;
    vector<int> first_child(dimension, -1);
    vector<int> second_child(dimension, -1);

    struct pair
    {
        int first_idx;
        int second_idx;
    };

    vector<pair> map;

    do {
        k1 = rand() % (dimension -1);
        k2 = rand() % (dimension -1);  //in case of hitting last index, while loop won't break (numbers from 1 to n-1)
    } while (k1 == k2);

    if (k1 > k2) {
        swap(k1, k2);
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
        for (int j = 0; j < map.size(); j++) {
            if (first_parent[i] == map[j].first_idx || first_parent[i] == map[j].second_idx) {
                exists_in_map = true;
                break;
            }
        }

        if (first_child[i] == -1
            && first_child.end () == find(first_child.begin(), first_child.end(), first_parent[i])
            && !exists_in_map) {
            first_child[i] = first_parent[i];
            cout << endl << "I: " << i << " val: " << first_parent[i] << endl;
        }
    }

    for (int i = 0; i < dimension; i++) {
        int val = first_parent[i];
        if (first_child[i] == -1) {
            vector<int> used_values;
            used_values.push_back(val);
            while (first_child.end() != find(first_child.begin(), first_child.end(), val)) {
                for (int j = 0; j < map.size(); j++) {
                    if (val == map[j].first_idx &&
                        used_values.end() == find(used_values.begin(), used_values.end(), map[j].second_idx)) {
                        val = map[j].second_idx;
                        used_values.push_back(val);
                        break;
                    } else if (val == map[j].second_idx &&
                               used_values.end() == find(used_values.begin(), used_values.end(), map[j].first_idx)) {
                        val = map[j].first_idx;
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
        for (int j = 0; j < map.size(); j++) {
            if (second_parent[i] == map[j].first_idx || second_parent[i] == map[j].second_idx) {
                exists_in_map = true;
                break;
            }
        }

        if (second_child[i] == -1
            && second_child.end () == find(second_child.begin(), second_child.end(), second_parent[i])
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
                for (int j = 0; j < map.size(); j++) {
                    if (val == map[j].first_idx  &&
                        used_values.end() == find(used_values.begin(), used_values.end(), map[j].second_idx)) {
                        val = map[j].second_idx;
                        used_values.push_back(val);
                        break;
                    } else if (val == map[j].second_idx  &&
                               used_values.end() == find(used_values.begin(), used_values.end(), map[j].first_idx)) {
                        val = map[j].first_idx;
                        used_values.push_back(val);
                        break;
                    }
                }
            }
            second_child[i] = val;
        }
    }
}

void Graph::update_pheromones(vector<vector<double>> &pheromones, vector<vector<int>> &routes) {
    double q = dimension; // pheromone amount left on the route
    double ro = 0.5;

    for (int i = 0; i < routes.size(); i++) { // for each route repeat
        int route_for_i = calculate_cost(&routes[i][0]); // calculate ant i's route cost
        for (int j = 0; j < routes.size() - 1; j++) {
            int city = routes[i][j]; // j city from i ant
            int next_city = routes[i][j + 1];

            // updating pheromone values on the edges between two cities
            pheromones[city][next_city] = (1 - ro) * pheromones[city][next_city] + q / (double) route_for_i;
            pheromones[next_city][city] = (1 - ro) * pheromones[next_city][city] + q / (double) route_for_i;
        }
    }

}

double Graph::get_phi(int first_city, int second_city, Ant *ant, vector<vector<double>> &pheromones) {
    double a = 1.1;
    double beta = 5.5; // beta parameter adjusting the influence of visibility

    // eta of going to another city
    double eta_ij = (double) pow(1.0 / matrix[first_city][second_city], beta);
    double tau_ij = (double) pow(pheromones[first_city][second_city], a);
    double sum = 0;

    for (int i = 0; i < dimension; ++i) {
        if (i == first_city)
            continue;
        if (!ant->visited[i]) {
            double eta = (double) pow(1.0 / matrix[first_city][i], beta);
            double tau = (double) pow(pheromones[first_city][i], a);
            sum += eta * tau;
        }
    }

    return (eta_ij * tau_ij) / (sum);
}

int Graph::get_next_city(vector<double> &probabilities) {
    double val = (double) rand() / ((double) (RAND_MAX + 1));

    int i = 0;
    double sum = probabilities[i];
    while (sum < val) { // choose next city where the sum of probabilities meets the random value
        i++;
        sum += probabilities[i];
    }

    return i;
}

void Graph::calculate_ant_routes(Ant *ant, vector<vector<int>> &routes, vector<vector<double>> &pheromones) {
    vector<double> probabilities;

    routes[ant->number][0] = ant->number;   // wierzcholek poczatkowy (0) dla mroki n to n
    ant->visited[ant->number] = true;       // jest odwiedzony

    for (int i = 0; i < dimension - 1; ++i) {
        int city_i = routes[ant->number][i];    //i-te miasto mrowki
        probabilities.clear();
        probabilities.resize(dimension, 0.0);
        for (int second_city = 0;
             second_city < dimension; second_city++) {    // liczenie prawdopodobienstwa dla wszsystkich miast
            if (city_i == second_city)
                continue;
            if (!ant->visited[second_city]) {
                probabilities[second_city] = get_phi(city_i, second_city, ant, pheromones);
            }
        }
        routes[ant->number][i + 1] = get_next_city(probabilities);    // zdecydowanie do ktorej krawedzi pojsc
        ant->visited[routes[ant->number][i + 1]] = true;  // zaznaczenie w liscie odwiedzonych
    }
}

void Graph::ant_colony_optimization() {
    int result = INT_MAX;

    int iteration_number = 100;
    int number_of_ants = dimension;
    vector<vector<int>> ant_routes(number_of_ants);
    vector<vector<double>> pheromones(dimension);

    chrono::steady_clock::time_point begin = chrono::steady_clock::now();

    for (int i = 0; i < number_of_ants; ++i) {
        ant_routes[i].resize(number_of_ants, -1);
    }

    for (int i = 0; i < dimension; ++i) {
        pheromones[i].resize(dimension);
        for (int j = 0; j < dimension; ++j) {
            pheromones[i][j] = (double) rand() / (double) RAND_MAX * dimension / matrix[0][1];
        }
    }

    for (int i = 0; i < iteration_number; ++i) {
        for (int j = 0; j < number_of_ants; ++j) {
            for (vector<int>::iterator it = ant_routes[j].begin(); it != ant_routes[j].end(); it++) {
                *it = -1; // przygotowanie trasy
            }
            Ant *ant = new Ant(j, dimension);
            calculate_ant_routes(ant, ant_routes, pheromones);
        }
        update_pheromones(pheromones, ant_routes);
    }


    for (int k = 0; k < dimension; ++k) {
        int temp = calculate_cost(&ant_routes[k][0]);
        if (temp < result) {
            result = temp;
            best_path = &ant_routes[k][0];
        }
    }

    cout << "Best Ant Colony Optimization result: " << result << "\n";
    cout << "Best path: ";
    for (int i = 0; i < dimension; i++) {
        if (i < dimension - 1) {
            cout << best_path[i] << ", ";
        }
        else {
            cout << best_path[i] << "\n";
        }
    }
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    cout << "Duration time: " << chrono::duration_cast<chrono::milliseconds>(end - begin).count() << "[ms]" << "\n";
    cout << "\n";

}

void Graph::ga_choose_params(move_foo *move, neighbourhood_type *nt, move_linear_update *linear_update,
        crossover_type *crossover) {
    char option;

    cout << "\n";
    cout << "Choose neighbourhood type (swap by default): \n";
    cout << " 1 - Swap.\n";
    cout << " 2 - Invert.\n";
    cout << " 3 - Insert.\n";
    cout << "\n";
    cin >> option;

    switch (option) {
        case '1':
            *move = swap;
            *linear_update = swap_linear_update;
            *nt = SWAP;
            break;
        case '2':
            *move = invert;
            *linear_update = invert_linear_update;
            *nt = INVERT;
            break;
        case '3':
            *move = insert;
            *linear_update = insert_linear_update;
            *nt = INSERT;
            break;
        default:
            *move = swap;
            break;
    }

    cout << "\n";
    cout << "Choose crossover type (OX by default): \n";
    cout << " 1 - OX.\n";
    cout << " 2 - PMX.\n";
    cout << "\n";
    cin >> option;

    switch (option) {
        case '1':
            *crossover = OX;
            break;
        case '2':
            *crossover = PMX;
            break;
        default:
            *crossover = OX;
            break;
    }
}


void Graph::genetic_alrogithm() {
    generate_population();
    int result = INT32_MAX;
    int first_rand_idx;
    int second_rand_idx;
    move_foo move; // mutation
    move_linear_update linear_update;
    neighbourhood_type ngbh_type;
    crossover_type crossover;

    ga_choose_params(&move, &ngbh_type, &linear_update, &crossover);
    chrono::steady_clock::time_point begin = chrono::steady_clock::now();

    best_path = &population[0][0];
    for (int j = 0; j < NUMBER_OF_GENERATIONS; ++j) {
        select_mating_pool_tournament();
        for (int i = 0; i < (int) (POPULATION_SIZE * CROSS_RATE); ++i) {
            int first_rand_idx, second_rand_idx;
            do {
                first_rand_idx = rand() % POPULATION_SIZE;
                second_rand_idx = rand() % POPULATION_SIZE;
            } while (first_rand_idx == second_rand_idx);

            crossover(population[first_rand_idx], population[second_rand_idx]);
        }

        for (int i = 0; i < (int) (POPULATION_SIZE * MUTATION_RATE); ++i) {
            int rand_index = rand() % POPULATION_SIZE;
            do {
                first_rand_idx = rand() % dimension;
                second_rand_idx = rand() % dimension;
            } while (first_rand_idx == second_rand_idx);
            move(&population[rand_index][0], first_rand_idx, second_rand_idx);
        }

        for (int k = 0; k < POPULATION_SIZE; ++k) {
            fitness[k] = calculate_cost(&population[k][0]);
            if (result > fitness[k]) {
                result = fitness[k];
                best_path = &population[k][0];
            }
        }
    }

    cout << "Best Genetic Algorithm result: " << result << "\n";
    cout << "Best path: ";
    for (int i = 0; i < dimension; i++) {
        if (i < dimension - 1) {
            cout << best_path[i] << ", ";
        }
        else {
            cout << best_path[i] << "\n";
        }
    }
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    cout << "Duration time: " << chrono::duration_cast<chrono::milliseconds>(end - begin).count() << "[ms]" << "\n";
    cout << "\n";
}