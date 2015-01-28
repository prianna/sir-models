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


int main(int argc, const char * argv[])
{
    for (int i = 0; i < MAX_RUNS; i++)
    {
        CrossScaleSim sim;
        
        // Note: the program will attempt to loop for number of generations
        // specified, but if transmission chain goes extinct prior to reaching max
        // generation, program will continue to run. Need to fix this.
        do
        {
            sim.cycle();
        }
        while ( sim.get_cycle() < MAX_CYCLES );
        //sim.calculateR();
        sim.output_map();
        //std::cout << "R_0 for this run is: " << sim.get_R0() << std::endl;
        sim.output_formatted();
        //sim.output_data();
    }
   
    
}





