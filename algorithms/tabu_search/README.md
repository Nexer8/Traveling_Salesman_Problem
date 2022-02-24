# Tabu Search

Tabu Search is a metaheuristic used to solve optimization problems. It is used in the search for optimal solutions or not much different from them for problems in various domains. It is based on an iterative search of the space of a solution for a dynamically changing neighborhood of a given solution and remembering some movements and their frequency in order to avoid getting stuck in local minima and looking for globally optimal solutions in an acceptable time.

## Implemented elements specific to the Tabu Search method

- **Number of iterations**: the number of iterations is a parameter that can be changed before executing the algorithm.
- **Tabu list length / cadency**: cadency is the number of iterations that a specific movement remains on the tabu list - is banned - and is also a parameter that is specified before the algorithm is run.
- **Different ways of determining the initial solution**
  - The cost of the natural solution, i.e. one for which the cities in the path are in the same order as in the distance matrix, is calculated as using the `calculateCost()` function, which takes as an argument a reference to a vector containing the path (the order of the cities that the salesman visits).
  - The `geRandomPathCost()` function takes a reference to a vector with path as an argument, and then randomizes it. The return value of this function is the cost of the newly generated path.
  - The `greedy()` function returns the cost of the path generated using the *greedy heuristic*. As arguments, the function takes a refence to the table with the path and information about the starting point, which is the first city the traveling salesman starts from.
- **Different neighborhood operations**
  - The **swap** operation changes the order of two components in the path.
  - The **inversion** operation involves reversing the order of the components on the route between the indices supplied as function parameters.
  - The **insert** operation inserts an element from the vector representing the path at the index specified in the second parameter into the location specified in the first argument. As a result, depending on the value of the parameters, the vector expands or contracts.
- **Tabu list representation**: it has been implemented as a vector of vectors, i.e. a matrix corresponding to the distance matrix taken from the file. So it has a constant size *NxN*, where *N* is the number of cities the salesman must pass through before returning to the starting city.
- **The aspiration criterion**: if the solution in a given neighborhood is better than the current best global solution, then the best global solution is updated without checking if the traffic is already present in tabu list, which is also updated - the reverse move will be forbidden for the number of iterations specified by the *cadency* term.
- **Critical event**: defined as the number of iterations after which a new starter solution should be generated and the tabu list cleared and the tabu search algorithm should be rerun.
- **Enhanced speed in neighborhood operations**: it can be noticed that counting the cost of a given path over time *O(n)* is not necessary because it is known how it will be changed. Therefore, it is possible to update the cost in linear time *O(1)* by removing the changed edges and adding new ones in their place. These enhancements are to be found inside `swap_linear_update()`, `invert_linear_update()` and `insert_linear_update()` functions.
