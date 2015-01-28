//
//  WithinHostODE.h
//  CrossScaleSim
//
//  Created by Prianna Ahsan on 7/11/13.
//  Copyright (c) 2013 Prianna Ahsan. All rights reserved.
//

#ifndef __CrossScaleSim__WithinHostODE__
#define __CrossScaleSim__WithinHostODE__

#include <vector>
#include "globals.h"
#include "odeint.hpp"


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
    /*
    // Initialise the equations and Runge-Kutta integration
    void Diff(double Pop[K]);
                                
    void Runge_Kutta(double step);
    //void Runge_Kutta_Fehlberg( double step);
    */
    std::vector<double> V, dVdt; // Vector of viral populations.
    
    double t, t_final;
    double step;


};

#endif /* defined(__CrossScaleSim__WithinHostODE__) */
