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
#include "ODESolver.h"

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
    
    std::vector<double> get_R_typeK()
    {
        return R_typeK;
    }
    
    std::vector<double> get_R_mixedK()
    {
        return R_mixedK;
    }
    
    double get_R_both()
    {
        return R_both;
    }
    
private:
    // Private classes used to hold dR0, integrate over this to compute R0 for a single strain.
    // Should only be used by CalculateR.
    class R0
    {
    public:
        R0( int i )
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
    
    class R0Both
    {
    public:
        R0Both( double T )
        :t_0(T)
        {
            for( int i = 0 ; i < K ; ++i )
            {
                V.push_back(V_CRIT);
            }
            
            system.set_V(V);
        }
        
        double operator()( double t ) const
        {
            V = system.IterateR( t_0, t );
            std::vector<double> integrand;
            integrand.push_back(1.0);
            for( int i = 0; i < K; ++i)
            {
                double expR = 1.0/exp( B[i]*V[i] );
                expR = 1.0-expR;
                integrand[0] *= expR;

            }
            return integrand[0];
        }
            
        private:
            double t_0;
            mutable std::vector<double> V;
            ODESolver system;
            
        };
    
    class R0TypeK
    {
    public:
        R0TypeK( int k )
        :k(k)
        {
            for( int i = 0 ; i < K ; ++i )
            {
                V.push_back(V_CRIT);
            }
            
            system.set_V(V);
        }
        
        double operator()( double t ) const
        {
            V = system.IterateR( T_0, t );
            std::vector<double> integrand;
            integrand.push_back( 1-( 1.0/exp( B[k]*V[k] )));
            for( int i = 0; i < K; ++i)
            {
                if( i != k )
                {
                    double expR = 1.0/exp( B[i]*V[i] );
                    integrand[0] *= expR;
                }
            }
            return integrand[0];
        }
        
    private:
        int k;
        mutable std::vector<double> V;
        ODESolver system;

    };
    
    class R0MixedK
    {
    public:
        R0MixedK( int k )
        :k(k)
        {
            for( int i = 0 ; i < K ; ++i )
            {
                V.push_back(V_CRIT);
            }
            
            system.set_V(V);
        }
        
        double operator()( double t ) const
        {
            V = system.IterateR( T_0, t );
            double integrand = 1 - ( 1.0/exp( B[k]*V[k] ));
            return integrand;
        }
        
    private:
        int k;
        mutable std::vector<double> V;
        ODESolver system;
        
    };

    
    // Calculates R0 by integrating Rfunction (below) from 0 to T.
    void calculateR();
    
    // Approximates R0 using Claude's approximation.
    void approxR();
    
    void calculateRMixedK();
    
    void calculateRBoth();
    
    void calculateRTypeK();
    
    std::vector<double> R_0, R_approx, R_typeK, R_mixedK;
    double R_both;
    
};



#endif
