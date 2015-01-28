//
//  globals.h
//  HybridSim
//
//  Created by Prianna Ahsan on 5/27/13.
//  Copyright (c) 2013 Prianna Ahsan. All rights reserved.
//

#ifndef HybridSim_globals_h
#define HybridSim_globals_h

const int MAX_RUNS = 1; // How many times the simulation will be run.
const double MAX_TIME=1; // How long each simulation will run for.

const double beta1=0.5; // Transmission term for strain 1.
const double beta2 = 0.6; // Transmission term for strain 2.
const double mu=0.0005; // Mutation rate for both strains.
const double lambda=9.1; // Constant birth rate.
const double dr = 0.05; // Constant death rate.
const double del=0.1; // Death rate of infected cells.
const double S0=50000.0; // Uninfected cells.
const double I10=10.0; // Cells infected with strain 1.
const double I20=10.0; // Cells infected with strain 2.

const int NUM_EVENTS = 5; // Number of events in sample space.

const int POP_THRESHOLD = 1000; // Infected pop threshold. When reached, program
                                // switches to purely deterministic simulation.

#endif
