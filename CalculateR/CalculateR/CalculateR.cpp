//
//  CalculateR.cpp
//  2StrainCSS
//
//  Created by Prianna Ahsan on 9/23/13.
//  Copyright (c) 2013 Prianna Ahsan. All rights reserved.
//

#include "CalculateR.h"

CalculateR::CalculateR()
{
    calculateR();
    approxR();
}

void CalculateR::calculateR()
{
    for (int i = 0; i < K; i++)
    {
        Rfunction dR(i);
        double R_temp = 1 - ( DEIntegrator<Rfunction>::Integrate(dR, T_0, T, 1e-5) )/T;
        R_temp = BETA_0*T*(R_temp);
        R_0.push_back( R_temp );
    }
}

void CalculateR::approxR()
{
    for (int i = 0; i < K; i++)
    {
        R_approx.push_back( BETA_0*( T + log( B[i]*V_CRIT )/RATES[i] ) );
    }
}