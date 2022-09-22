Parallel NSGAII
========================================================================
This code has been modified from the original NSGA-II code 
developed by Kalyanmoy Deb so that it can be executed in parallel
on HPC systems for costly problem evaluations. It uses OpenMP for fine scale parallelism and MPI for coarse scale parallelism.
Additionally, a user can specify an external shared library to be used for objective functions and constraints 
calculations.


About the Algorithm
--------------------------------------------------------------------------
NSGA-II: Non-dominated Sorting Genetic Algorithm - II

Please refer to the following paper for details about the algorithm:

Authors: Dr. Kalyanmoy Deb, Sameer Agrawal, Amrit Pratap, T Meyarivan
Paper Title: A Fast and Elitist multi-objective Genetic Algorithm: NSGA-II
Journal: IEEE Transactions on Evolutionary Computation (IEEE-TEC)
Year: 2002
Volume: 6
Number: 2
Pages: 182-197

The original source code was developed by [Dr. Kalyanmoy Deb](http://www.iitk.ac.in/kangal)

How to compile and run the program
---------------------------------------------------------------------------
Makefile has been provided for compiling the program on linux (and unix-like)
systems. Edit the Makefile to suit your need. By default, provided Makefile
attempts to compile and link all the existing source files into one single
executable. Comment out the USE_MPI option to disable MPI. Comment out the USE_OPENMP option to disable OpenMP.

The command to use for compilation under Ubuntu 18 is:

	make -f Makefile.chpc

Name of the executable produced is: parallelnsga2r

To run the program type the following:
    
	 mpirun -n x ./parallelnsga2r random_seed <inp_file.in> -p


where x is the number of MPI processes to use and random_seed is a real number in (0,1) which is used as a seed for random number generator. <inp_file.in> is the problem input file. -p is an optional argument that indicates whether you want to write each individual in the population for each generation to file. Writing for each generation can be costly in terms of simulation times.

About the input parameters
---------------------------------------------------------------------------
1. popsize: This variable stores the population size (a multiple of 4)
2. ngen: Number of generations
3. nobj: Number of objectives
4. ncon: Number of constraints
5. nreal: Number of real design variables
6. lines (6->6+nreal)min_realvar[i]: minimum value of 'i' real variable  max_realvar[i]: maximum value of 'i' real variable (space separation)
7. pcross_real: probability of crossover of real variable
8. pmut_real: probability of mutation of real variable It is recommended each decision variable is mutated with a probability of 1 / L, where L is the number of decision variables. This results in one mutation per offspring on average.
9. eta_c: distribution index for real variable SBX crossover
10. eta_m: distribution index for real variable polynomial mutation
11. nbin: number of binary variables
12. nbits[i]: number of bits for i^{th} binary variable
13. (lines 13+nreal->13+nreal+nbits) min_binvar[i]: minimum value of i^{th} binary variable max_binvar[i]: maximum value of i^{th} binary variable (space separation)
14. pcross_bin: probability of crossover for binary variable
15. pmut_bin: probability of mutation for binary variable
16. problem definition: -t x for one of the test problems or -f <shared_library.in> funcname for a custom problem. Where x in is the index for the test problem and  <shared_library.in> is the path to the shared library, and funcname is the name of the function for the problem definition. The following are the indexes for the test problems
    * 0 = SCH1
    * 1 = SCH2
    * 2 = FON
    * 3 = KUR
    * 4 = POL
    * 5 = VNT
    * 6 = ZDT1
    * 7 = ZDT2
    * 8 = ZDT3
    * 9 = ZDT4
    * 10 = ZDT5
    * 11 = ZDT6
    * 12 = BNH
    * 13 = OSY
    * 14 = SRN
    * 15 = TNK
    * 16 = CTP1
    * 17 = CTP2
    * 18 = CTP3
    * 19 = CTP4
    * 20 = CTP5
    * 21 = CTP6
    * 22 = CTP7
    * 23 = CTP8
    * 24 = SCHF : SIX-HUMP CAMEL FUNCTION (https://www.sfu.ca/~ssurjano/camel6.html) 
        sol. f(x^*)=-1.0316 at x^*=(0.0898,-0.7126) \& (-0.0898, 0.7126)
17. remaining lines: will be read as an array of strings (vector<string>&) and passed as an argument to the problem definition

About the output files
---------------------------------------------------------------------------
* initial_pop.out: This file contains all the information about initial population.
* final_pop.out: This file contains the data of final population.
* all_pop.out: This file containts the data of populations at all generations.
* best_pop.out: This file contains the best solutions obtained at the end of simulation run.
* params.out: This file contains the information about input parameters as read by the program.

How to create your own problem definition
---------------------------------------------------------------------------
The files "problemdef.h & problemdef.c" contain the problem header and definition. First copy one of the definitions of the header and change its name, define it at the begining of the .c file and create a new function defintion. The variables used in these functions are:
	
* xreal[i]		: 
* nobj        		: number of objective values
* obj[i]      		: 'i' objective value
* ncon        		: number of constraint values
* constr[i]   	: 'i' constraint value
* optionalArgs	: Other arguments passed by input file in vectorial string format (vector<string>&)	
	

The following points need to be kept in mind while writing the objective and constraint
functions.
1. The code has been written for minimization of objectives (min f_i). If you want to
maximize a function, you may use negetive of the function value as the objective value.
2. A solution is said to be feasible if it does not violate any of the constraints.
Constraint functions should evaluate to a quantity greater than or equal to zero
(g_j >= 0), if the solution has to be feasible. A negative value of constraint means that the constraint
it is being violated, it is not a feasible solution.
3. If there are more than one constraints, it is advisable (though not mandatory)
to normalize the constraint values by either reformulating them or dividing them
by a positive non-zero constant (example TODO). 

Some sample problems (24 test problems from Dr. Deb's book - Multi-Objective Optimization
using Evolutionary Algorithms) have been provided as examples to guide you
define your own objective and constraint functions. You can also link other
source files with the code depending on your need.

About the files
---------------------------------------------------------------------------
    global.h: Header file containing declaration of global variables and functions
    rand.h: Header file containing declaration of variables and functions for random number generator
    allocate.c: Memory allocation and deallocation routines
    auxiliary.c: auxiliary routines (not part of the algorithm)
    crossover.c: Routines for real and binary crossover
    crowddist.c: Crowding distance assignment routines
    decode.c: Routine to decode binary variables
    dominance.c: Routine to perofrm non-domination checking
    eval.cpp: Routine to evaluate constraint violation
    fillnds.c: Non-dominated sorting based selection
    initialize.c: Routine to perform random initialization to population members
    list.c: A custom doubly linked list implementation
    merge.c: Routine to merge two population into one larger population
    mutation.c: Routines for real and binary mutation
    parallelnsga2r.cpp: Implementation of main function and the NSGA-II framework
    problemdef.c: Test problem definitions
    rand.c: Random number generator related routines
    rank.c: Rank assignment routines
    report.c: Routine to write the population information in a file
    sort.c: Randomized quick sort implementation
    tourselect.c: Tournament selection routine



---------------------------------------------------------------------------
Contact [me](caleb.buahin@usu.edu) with questions or comments.


[highlight.js]: http://softwaremaniacs.org/soft/highlight/en/