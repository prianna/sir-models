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
const int MAX_RUNS = 10000;

// Size of viral type vector
const int K = 1;
//const int K = 2;


// Starting viral pop
const std::vector<int> VPOP = {100};
//const std::vector<int> VPOP = {100, 0};

// Growth rate for within-host phase.
const double RATES[K] = {1.0};
//const double RATES[K] = {1.0, 1.0};

// Rel. transmissibility
const double B1 = pow(10.0, 10.0);
//const double B2 = pow(10.0, 9.0);
//const double B[K] = {1/B1, 1/B2};
const double B[K] = {1/B1};

// Length of infection period after viral load hits V_CRIT.
const double T = 15.0;
const double T_0 = 0.0;
const double T_max = 15.0;


// Contact rate
const double BETA_0 = 19.4116/T;

// Mutation rate
const double denom = pow(10.0, 8.0);
const double mu = 0.0;
//const double mu = 1.0/denom;

// Max number of generations to run simulation.
const int MAX_GENS = 1;

// Max overall number of individuals allowed in the population.
const int POP_MAX = 1000;

// Critical threshold for infection to occur.
const double V_CRIT = pow(10.0, 4.0);

// Scaling factor. If particles exceed N, we scale back by N.
const double N = pow(10.0, 10.0);


#endif