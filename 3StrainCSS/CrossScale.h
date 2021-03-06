//
//  CrossScale.h
//  CrossScaleSim class definition.


#ifndef __CrossScaleSim__CrossScale__
#define __CrossScaleSim__CrossScale__

#include "globals.h"
#include "ODESolver.h"
#include <random>
#include <chrono>
#include <cmath>
#include <map>
#include <list>
#include <iostream>
#include <iomanip>
#include <vector>
#include <sstream>
#include <fstream>
#include <queue>


class CrossScaleSim
{
public:
    // Constructor.
    CrossScaleSim();
    
    // Destructor.
    ~CrossScaleSim()
    {
        delete population.front();
    }
    
    // Main iterative loop.
    void cycle()
    {
        if (population.size() >= POP_MAX)
        {
            cur_cycle = MAX_CYCLES;
            return;
        }
    
        else if ( !children.empty() )
        {
            while ( !children.empty() && generations <= MAX_GENS )
            {
                individual *parent = children.front();
                generations = parent->gen;
                double v_load = critical(parent);
                if ( v_load < V_CRIT )
                {
                    std::vector<int> v_crit;
                    ODESolver initial(parent->v_0, T_0, T);
                    v_crit = initial.Iterate( 1.0 );
                    parent->t_crit = parent->t_0+initial.get_time();
                    transmission(v_crit, parent);
                }
                else
                    transmission(parent->v_0, parent);
                
                population.push_back(parent);
                children.pop();
            }
            cur_cycle = MAX_CYCLES;
        }
        cur_cycle++;
    }
    
    
    // Accessor.
    int get_cycle()
    {
        return cur_cycle;
    }
    
    // Calculation of R_naught
    void calculateR();

    double get_R0()
    {
        return R_0;
    }
    
    // Outputs a table of individuals and all relevant data to the console.
    void output_map();
    
    // Outputs the following: time created (t_0), viral load at t_0, generation
    // parent to respective files. 
    void output_data();
    
    void output_formatted();

private:
    
    // Node to contain data about individuals. Private struct, only accessed by
    // members of class CrossScaleSim.
    struct individual
    {
        double t_0; // Time of creation.
        double t_crit; // Time to v.crit.
        std::vector<int> v_0; // Viral types present at t_0.
        individual *parent; // Pointer to parent.
        std::vector<individual*> children; // Vector of pointers to children.
        int key; // Unique identifier.
        int gen; // Generation created.

        
        // default constructor, takes 3 parameters: time of infection, viral types
        // at infection, pointer to parent.
        individual( double tc, std::vector<int> vir_type, individual *parent, int ID, int gen)
        : t_0(tc), v_0(vir_type), parent(parent), key(ID), gen(gen), t_crit(tc)
        {}
        
        // Destructor, to free up memory.
        ~individual()
        {
            for (int i = 0; i < children.size(); i++)
            {
                delete children[i];
            }
            children.clear();

        }
    };
    
    int critical( individual *parent )
    {
        int viral_load = 0.0;
        
        for (int i = 0; i < K; i++)
        {
            viral_load += parent->v_0[i];
        }
        
        return viral_load;
    }
        
    // Precondition: Time is a number drawn from a uniform distribution (0.0, T)
    // and vtypes is a vector containing the initial viral load of the parent.
    // Postcondition: A vector containing the viral load at t = time is returned.
    std::vector<double> withinHost( double time_init, double time_fin, std::vector<int> vtypes );

    // Precondition: Parentvir contains the initial viral load of the parent.
    // Postcondition: The individual undergoes its transmission cycle, and
    // infects c children, where c is generated by poisson(B0*T). Each child
    // is added to the individual's children vector.
    void transmission( std::vector<int> parentvir, individual* parent );
    
    std::uniform_real_distribution<double> unif; // Uniform distribution.
    std::mt19937 gen; // Mersenne Twister RNG.
    std::poisson_distribution<> poiss; // Poisson distribution.
    std::vector<individual*> population;
    std::queue<individual*> children;
    
    // Initial vector of viral loads.
    std::vector<int> V_0;
    int generations, cur_cycle;
    double R_0;
};

#endif /* defined(__CrossScaleSim__CrossScale__) */
