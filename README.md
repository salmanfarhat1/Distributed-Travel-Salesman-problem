# Distributed-Travel-Salesman-problem
Traveling Salesman Problem (TSP) is one of the most common studied problems in combinatorial optimization. Given the list of cities and distances between them, the problem is to find the shortest tour possible which visits all the cities in list exactly once and ends in the city where it starts.

Procedure that has been followed: • Arguments parsing • Reading input file • Population initialization • Fitness Calculation • Loop • Selection • Crossover • Mutation • Evaluation • Replacement

there are two versions of crossover: Greedy Crossover V1 For a pair of parents i & j • Pick 1st city from parent i • Next choose nearest next city not causing cycle from i (or else j) (i.e. not already figuring so far in chromosome) • If nearest next city from both parents is causing a cycle, then a random city (not introducing a cycle) is chosen • Loop until all cities are chosen

Greedy Crossover V2 For a pair of parents • Rand-pick a city as start • Next chosen city is the k-th nearest one not causing cycle (i.e. not already figuring so far in chromosome) • If nearest next city from both parents is causing a cycle, then a random city (not introducing a cycle) is chosen • Loop until all cities are chosen
