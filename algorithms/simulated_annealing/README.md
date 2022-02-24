# Simulated Annealing

Simulate Annealing is a metaheuristic approach based on the events that occur during the liquid cooling and metal cooling processes. This is a metallurgical process known as annealing. It is used to solve situations when finding a solution that is close to the optimal global solution is more essential than determining the ideal local solution accurately in a specific time.

## Implemented elements specific to the Simulated Annealing method

- **Setting the initial temperature**: the initial and final temperature, which are also the final condition of the algorithm, are parameters that can be changed before the algorithm is executed.
- **Various cooling schemes**: three different cooling methods (`linear_cooling()`, `geometric_cooling()`, and `logarithmic_cooling()`) have been implemented in the algorithm representing simulated annealing, i.e. methods determining the speed of reaching the final conditions, and thus the speed of the algorithm.
- **Different approaches to compute the initial solution**: different ways to compute the initial solution are implemented in the same way as in the *Tabu Search* algorithm.
- **Repeating iterations for a given temperature**: the value of the *repeat* variable determines the number of iterations that will be performed before the temperature changes.
- **Various neighborhood operations with enhanced speed**: see *Tabu Search* for more details.
- **Stagnation**: increases when the current solution is not better than the best solution found. When its value exceeds the set limit, the temperature value increases. In case of finding a value better than the best one so far, the stagnation indicator is reset.
- **Strength of the solution**: the probability of accepting a solution is expressed by the formula: $P(\Delta E) = \exp(\frac{-\Delta E}{T})$. It is implemented as the `get_probability()` function.
- **Deciding on the execution of the next step**: checking the condition whether the difference between the present and the best result found is less than 0 or, if it is not, checking whether the function returning the probability value gives values greater than those randomly generated in the range [0; 1]. If any of the above conditions are true, execute the given step.
