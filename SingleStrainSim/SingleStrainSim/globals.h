//
//  globals.h
//  Modify these to alter parameters and customize simulation.
//
//  Created by Prianna Ahsan on 6/19/13.
//  Copyright (c) 2013 Prianna Ahsan. All rights reserved.
//
#ifndef CrossScaleSim_globals_h
#define CrossScaleSim_globals_h
#include <cmath>
#include <vector>

// Max number of runs using parameters.
const int MAX_RUNS = 1000;

// Max number of attempts to try and generate transmission using current RNG.
const int MAX_CYCLES = 1;

// Size of viral type vector
const int K = 1;

// Starting viral pop
const std::vector<int> VPOP = {1};

// Growth rate for within-host phase.
const double RATES[K] = {1.0};

// Rel. transmissibility
const double B1 = pow(10.0, 6.3);
const double B[K] = {1/B1};

// Contact rate.
const double BETA_0 = 2.0/15;

// Mutation rate
const double denom = pow(10.0, 5.0);
const double mu = 1.0/denom;

// Length of infection period after viral load hits V_CRIT.
const double T = 15.0;
const double T_0 = 0.0;

// Max number of generations to run simulation.
const int MAX_GENS = 1;

// Max overall number of individuals allowed in the population.
const int POP_MAX = 10000;

// Critical threshold for infection to occur.
const double V_CRIT = 10000.0;

// Scaling factor. If particles exceed N, we scale back by N.
const double N = pow(10.0, 10.0);



#endif
