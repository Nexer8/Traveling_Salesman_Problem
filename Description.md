# TSP

The *travelling salesman problem* (also known as the *TSP*) asks, "Given a list of cities and the distances between each pair of cities, what is the shortest feasible route that visits each city precisely once and returns to the starting city?" It is an *NP-hard* problem, which means that it is difficult to solve for large inputs and there is no known algorithm that can solve it in polynomial time.

This repository contains various approaches to solving the *TSP* implemented in *C++*.

# Brute Force

The first one, and the easiest one is the brute force approach. It is a simple algorithm that tries all possible permutations of the cities and returns the shortest one. It is very slow, but it is a good starting point to understand the problem. It always yields the proper outcome. This, however, comes at the price of a computational complexity.

## Swap-based Brute Force

This approach used *swap* and recursion. Its computational complexity is *O(n!)*, resulting in an extremely long execution time.

What I did here is to swap the cities in the order of the path. The first city is fixed, and the rest are swapped. The algorithm is then recursively called for the rest of the cities.

## Tree search Brute Force

Another algorithm that can be used to solve small instances of the traveling salesman problem is the brute force based on the searched tree and recursion. The computational complexity in this case is the same as in the swap example and is *O(n!)*. The algorithm works slowly but guarantees the correct result.

# Branch and Bound

In the worst-case scenario, the computational complexity is identical to that of the brute force approach, i.e. *O(n!)*. In fact, however, this strategy works extremely rapidly, although it is dependent on the exact instance of the issue as well as the function used to determine upper and lower bounds. This particular implementation uses recursion.

# Dynamic Programming

Dynamic programming is a technique for addressing optimization problems that is based on the search for a decision sequence. This approach computes solutions to all sub-problems, beginning with the most basic. The outcomes of solved sub-problems are preserved in the constructed table, so in the event of running into the same issue again, the algorithm will not have to solve it again. The algorithm has an $O(n^22^n)$ time complexity and an $O(n2^n)$ memory complexity.

The Held-Karp algorithm is implemented here.

