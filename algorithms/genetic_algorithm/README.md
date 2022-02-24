# Genetic Algorithm

The Genetic Algorithm is a metaheuristic inspired by the process of natural selection. It is used to generate high quality solutions to optimization problems. The genetic algorithm belongs to the group of evolutionary algorithms. This metaheuristic uses operators inspired by the theory of evolution. These include selection, crossing and mutation.

## Implemented elements specific to the Genetic Algorithm

- **Population size**: a parameter that can be changed before executing the algorithm.
- **Number of generations**: the number of iterations the genetic algorithm is repeated.
- **Mating pool tournament**: the parameter determining the size of the tournament, i.e. the number of individuals taken into account in one tournament, can be changed before the algorithm is executed. The tournament selection itself consists in choosing the best individual out of all those taking part in the tournament.
- **Crossover operators**: three different crossover operators are implemented in the program. They serve to explore solved space, which means that the algorithm will be able to come out of the local minimum.
  - *Order Crossover (OX)*
  - *Partially-mapped Crossover (MPX)*
  - *Non-wrapping Ordder Crossover (NWOX)*
- **Mutation operators**: three different mutation operators are implemented in the program. The mutation operators serve to exploit the solution space. These operations are the same as the neighborhood operations (see *Tabu Search*).
- **Crossover probability**: the value of the probability of crossover is determined before starting the execution of the algorithm.
- **Mutation probability**: the value of the mutation probability is determined before starting the algorithm execution.
