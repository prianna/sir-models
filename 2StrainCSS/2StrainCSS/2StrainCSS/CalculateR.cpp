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
    calculateRBoth();
    calculateRTypeK();
    calculateRMixedK();
}

void CalculateR::calculateR()
{
    for (int i = 0; i < K; i++)
    {
        R0 dR(i);
        double R_temp = ( DEIntegrator<R0>::Integrate(dR, T_0, T, 1e-8) );
        R_temp = T-R_temp;
        R_0.push_back( ( BETA_0*R_temp )  );
    }
}


void CalculateR::calculateRBoth()
{
    R0Both dR(T_0);
    double R_temp = ( DEIntegrator<R0Both>::Integrate(dR, T_0, T, 1e-8) );
    R_both = R_temp*BETA_0;
    
}

void CalculateR::calculateRTypeK()
{
    for (int i = 0; i < K; i++)
    {
        R0TypeK dR(i);
        double R_temp = ( DEIntegrator<R0TypeK>::Integrate(dR, T_0, T, 1e-8) );
        R_typeK.push_back(R_temp*BETA_0);
    }
}

void CalculateR::calculateRMixedK()
{
    for (int i = 0; i < K; i++)
    {
        R0MixedK dR(i);
        double R_temp = ( DEIntegrator<R0MixedK>::Integrate(dR, T_0, T, 1e-8) );
        R_mixedK.push_back(R_temp*BETA_0);
    }
}


void CalculateR::approxR()
{
    for (int i = 0; i < K; i++)
    {
        R_approx.push_back( BETA_0*( T + log( B[i]*V_CRIT )/RATES[i] ) );
    }
}