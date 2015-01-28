//
//  CalculateR.h
//  2StrainCSS
//
//  Created by Prianna Ahsan on 9/23/13.
//  Copyright (c) 2013 Prianna Ahsan. All rights reserved.
//

#ifndef _StrainCSS_CalculateR_h
#define _StrainCSS_CalculateR_h

#include "globals.h"
#include "odeint.hpp"
#include "DEIntegrator.h"

class CalculateR
{
public:
    CalculateR();
    
    std::vector<double> get_R0()
    {
        return R_0;
    }
    
    std::vector<double> get_R_approx()
    {
        return R_approx;
    }
    
    
private:
    // Private class used to hold dR0, integrate over this to compute R0 for a single strain.
    // Should only be used by CalculateR.
    class Rfunction
    {
    public:
        Rfunction( int i )
        :i(i)
        {
        }
        
        double operator()( double t ) const
        {
            double expR = exp( RATES[i]*t );
            double minim = fmin( V_CRIT*expR, N );
            expR = 1.0/(exp(B[i]*minim));
            return expR;
        }
        
    private:
        int i;
        
    };
    
    // Calculates R0 by integrating Rfunction (below) from 0 to T.
    void calculateR();
    
    // Approximates R0 using Claude's approximation.
    void approxR();
    
    std::vector<double> R_0, R_approx;
    
    
};



#endif
