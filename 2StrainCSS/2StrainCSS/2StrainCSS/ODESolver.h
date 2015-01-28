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
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

typedef std::vector< double > state_type;

class ODESolver
{
public:
    ODESolver( std::vector<int> init_types, double t_init, double t_fin );
    
    ODESolver();
    
    std::vector<double> Iterate();
    
    std::vector<int> Iterate( double time );
    
    std::vector<double> IterateR( double T0, double TF) const;
    
    double get_time()
    {
        return t;
    }
    
    void set_V( std::vector<double> types)
    {
        for (int i = 0; i < K; i++)
        {
            V.push_back(types[i]);
        }
    }

private:
    /*
    // Initialise the equations and Runge-Kutta integration
    void Diff(double Pop[K]);
                                
    void Runge_Kutta(double step);
    //void Runge_Kutta_Fehlberg( double step);
    */
    void prepare_output( std::ofstream &output, std::string filename );
    void end_output( std::ofstream &output )
    {
        output.close();
    }
    
    struct streaming_observer
    {
        std::ofstream &m_out;
        streaming_observer( std::ofstream &out ) : m_out( out ) {}
        
        void operator()( const state_type &x , double t ) const
        {
            m_out << t;
            for( size_t i=0 ; i < x.size() ; ++i )
            {
                m_out << ", " << x[i];
            }
            m_out << "\n";
        }
    };

    std::vector<double> V, dVdt; // Vector of viral populations.
    std::ofstream formatted_output;
    double t, t_final;
    double step;
    


};

#endif /* defined(__CrossScaleSim__WithinHostODE__) */
