//
//  simulations.h
//  HybridSim
//
//  Created by Prianna Ahsan on 5/27/13.
//  Copyright (c) 2013 Prianna Ahsan. All rights reserved.
//

#ifndef __HybridSim__simulations__
#define __HybridSim__simulations__

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <cstring>
#include <sstream>
#include <cmath>
#include <ctime>
#include "globals.h"

class Simulation
{
public:
    Simulation()
    : events(0), t(0.0), I1(I10), I2(I20), S(S0), dt(0)
    {
    }
    
    // Accessor functions
    double getS(){ return S; }
    double getI1(){ return I1; }
    double getI2(){ return I2; }
    double getdt(){ return dt; }
    double gett(){ return t; }
    
    // Mutator function
    void Iterate( FILE *output );

private:
    
    // Runs if either infected population is < 1000.
    void Runge_Kutta( double timeStep );
    
    // Runs if either infected population size is > 1000.
    void Runge_Kutta( double timeStep, int n );
    
    void Diff( double pop[3] );
    double Gillespie();
    
    void Output_Data( FILE *output );
    
    double t, S, I1, I2, dt;
    double dPop[3];
    int events;
};




#endif /* defined(__HybridSim__simulations__) */
