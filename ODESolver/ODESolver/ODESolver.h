//
//  ODESolver.h
//  ODESolver
//
//  Created by Prianna Ahsan on 8/13/13.
//  Copyright (c) 2013 Prianna Ahsan. All rights reserved.
//

#ifndef __ODESolver__ODESolver__
#define __ODESolver__ODESolver__

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <cmath>

//Globals
const int K = 2;
const double mu = 1/pow(10, 5);

class ODESolver
{
public:
    ODESolver( std::vector<int> init_types, double t_init, double t_fin );
    
    std::vector<double> Iterate();
    
    std::vector<int> Iterate( double totalLoad );
    
    double get_time()
    {
        return t;
    }
    
private:
    // Initialise the equations and Runge-Kutta integration
    void Diff(double Pop[K]);
    
    void Runge_Kutta(double step);
    
    void Runge_Kutta_Fehlberg( double step );
    
    double mu; // Mutation rate for both strains.
    std::vector<double> V; // Vector of viral populations.
    
    double t, t_final;
    double dPop[K];
    double step;
    
    
};


#endif /* defined(__ODESolver__ODESolver__) */
