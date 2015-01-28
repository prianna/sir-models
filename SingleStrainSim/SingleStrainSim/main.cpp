// Modify MAX_GENS in globals.h to change number of generations.
// See globals.h for other parameters.
// Possibile change: Add another function to allow for numerous cycles
// thus running variable separate simulations repeatedly with different
// starting conditions. Would move everything inside main into the new
// "inner" function. Each simulation would have a different random
// seed.


#include <iostream>
#include <string>
#include <sstream>
#include "CrossScale.h"
#include "CalculateR.h"


int main(int argc, const char * argv[])
{
    for (int i = 0; i < MAX_RUNS; i++)
    {
        CrossScaleSim sim_temp;
        
        // Note: the program will attempt to loop for number of generations
        // specified, but if transmission chain goes extinct prior to reaching max
        // generation, program will continue to run. Need to fix this.
        do
        {
            sim_temp.cycle();
        }
        while ( sim_temp.get_cycle() < MAX_CYCLES );
        
        
        sim_temp.output_map();
        sim_temp.output_formatted();
        int max_gen = sim_temp.get_gens();
        if ( max_gen >= MAX_GENS)
        {
            std::cout << "Max generation attained: " << max_gen << std::endl;
        }
        //sim.output_data();
    }
    CalculateR R;
    std::vector<double> R_0 = R.get_R0();
    std::vector<double> R_approx = R.get_R_approx();
    
    for (int i = 0; i < K; i++)
    {
        std::cout << "EXACT R_0_STRAIN" <<i<< " for these params is: " << R_0[i] << std::endl;
        std::cout << "APPROX R_0_STRAIN" <<i<< " for these params is: " << R_approx[i] << std::endl;

    }
}





