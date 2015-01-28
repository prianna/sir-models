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

// Max number of runs using parameters.
const int MAX_RUNS = 100;

// Max number of attempts to try and generate transmission using current RNG.
const int MAX_CYCLES = 10;

// Size of viral type vector
const int K = 3;

// Growth rate for exponential phase.
const double RATES[K] = {1.0, 1.0, 1.1};

const double B1 = pow(10.0, 10.0);
const double B2 = pow(10.0, 10.0);
const double B3 = pow(10.0, 6.0);

// Rel. transmissibility
const double B[K] = { 1/B1, 1/B2, 1/B3};

// Mutation rate
const double M1 = (pow(10, 5.0));
const double MU = (1/M1);

// Contact rate.
const double BETA_0 = 2.0/15;

// Max no of ticks per individual (time).
const double T = 15.0;

// Max number of generations to run simulation.
const int MAX_GENS = 10;

// Max overall number of individuals allowed in the population.
const int POP_MAX = 10000;

// Critical threshold for infection to occur.
const double V_CRIT = 10000.0;

// Scaling factor. If particles exceed N, we scale back by N.
const double N = pow(10.0, 10.0);

const double T_0 = 0.0;

const int V_INIT[K] = { 1, 100, 0};



#endif
