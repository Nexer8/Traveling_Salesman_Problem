# Ant Colony Optimization

The Ant Colony Optimization technique is a metaheuristic method that is based on the behavior of ants looking for food for their colony. The ants roam at random, and when they discover food, they return to their colony, leaving a trail of pheromones behind. It is a determinant, which means that when it encounters another ant, it begins to follow it in search of food. As a result, the ant stops traveling at random. However, the pheromones will dissipate with time. This allows to determine the best feeding path; otherwise, each succeeding route would have the same force as the preceding one. The algorithm's end condition might be a predetermined number of iterations or an occurrence in which all ants choose the same path.

## Implemented elements specific to the Ant Colony Optimization algorithm

- **Ant definition**: ant's attributes include its index and a vector in which the information whether a given city has already been visited by this ant is recorded.
- **Number of ants**: the number of ants was set to a value equal to the number of cities in the analyzed TSP problem.
- **Initial positions of ants**: the initial positions of the ants were set so that each ant would start its journey in a different city.
- **Initial setting of the number of pheromones on the edges**: initial setting the number of pheromones on the edges prevents ants from quickly selecting a *greedy* path.
- **Number of iterations**: the number of iterations through which the ant algorithm is repeated.
- **The way to set the ant's route**: the method of setting the ant routes is shown in the `calculateAntRoutes()` function.
- **Choosing the next city for the ant**: `getNextCity()` is the function returning the index of the newly selected city to which the ant will.
- **The probability of an ant going to a given city**: $p_{ij}=\begin{cases}{\frac{(\tau_{ij})^{\alpha}(\eta_{ij})^{\beta }}{\sum_{c_i, l \in \Omega}(\tau_{ic})^\alpha (\tau_{ic})^\beta}}\, ,& \bigvee c_{i,l} \in \Omega\\0\, ,& \bigvee c_{i,l} \notin \Omega \end{cases}$
  - *c* - another possible city,
  - *W* - acceptable solution,
  - $\eta_{ij}$ - value of the local visibility criterion function, where $\eta = \frac{1}{d{ij}}$,
  - *a* - parameter regulating the influence of $\tau_{ij}$,
  - *b* - parameter regulating the influence of $\eta_{ij}$.
- **Ways to update the amount of pheromones on the edge**: according to *Dorigo*
  - **Ant-Density (DAS)**: the amount of pheromone per length unit is always constant and independent of the edge length.
  - **Ant-Quantity (QAS)**: when crossing an edge, the constant amount of pheromone is divided by the length of this edge.
  - **Ant-Cycle (CAS)**: the fixed amount of pheromone is divided by the length of the route found by a given ant.
