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
#include <chrono>



int main(int argc, const char * argv[])
{    
   auto start = std::chrono::high_resolution_clock::now();
   
   for (int i = 0; i < MAX_RUNS; i++)
    {
        CrossScaleSim sim_temp;
        
        // Note: the program will attempt to loop for number of generations
        // specified, but if transmission chain goes extinct prior to reaching max
        // generation, program will continue to run. Need to fix this.
        sim_temp.cycle();
        
        sim_temp.output_map();
        sim_temp.output_formatted();
        //sim.output_data();
    }
    CalculateR R;
    std::vector<double> R_0 = R.get_R0();
    std::vector<double> R_approx = R.get_R_approx();
    std::vector<double> R_typeK = R.get_R_typeK();
    std::vector<double> R_mixedK = R.get_R_mixedK();
    double R_both = R.get_R_both();
    
    std::cout << std::endl << "RO (BOTH) " << " for these params is: " << R_both << std::endl;

    for (int i = 0; i < K; i++)
    {
        std::cout << "MIXED STRAIN R0 TYPE " <<i+1<< " for these params is: " << R_mixedK[i] << std::endl;
        std::cout << "2 STRAIN R0 TYPE " <<i+1<< " for these params is: " << R_typeK[i] << std::endl;
        std::cout << "EXACT R0 STRAIN " <<i+1<< " for these params is: " << R_0[i] << std::endl;
        std::cout << "APPROX R0 STRAIN " <<i+1<< " for these params is: " << R_approx[i] << std::endl;
        std::cout << B[i] << std::endl;
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Time elapsed: " << std::chrono::duration <double, std::milli> (elapsed).count() << " ms" << std::endl;
}





