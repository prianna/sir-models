//
//  CrossScale.cpp
//  Source code for CrossScaleSim class.

#include "CrossScale.h"


CrossScaleSim::CrossScaleSim()
: generations(0), cur_cycle(0), poiss(BETA_0*T), unif(0.0,T)
{
    std::random_device rd;
    gen.seed(rd()); // Seed Mersenne Twister pseduo RNG with random device.
	V_0 = VPOP;
    population.push_back(new individual(0.0, V_0, NULL, 1, generations));
    individual *parent = population.front();
    ODESolver initial(parent->v_0, T_0, T);
    std::vector<int> v_crit = initial.Iterate( 1.0 );
    parent->t_crit = initial.get_time();
    parent->t_abs = initial.get_time();
    transmission(v_crit, parent);
}

void CrossScaleSim::output_data()
{
    // Create ifstream objects to check if files exist.
    std::ifstream check_file_1 ("time.txt");
    std::ifstream check_file_2 ("virus.txt");
    std::ifstream check_file_3 ("gen.txt");
    std::ifstream check_file_4 ("parent.txt");
    std::ofstream file_1, file_2, file_3, file_4;

    // Check if file exists...
    if ( !check_file_1.is_open() ) // If the file doesn't exist, create it.
    {                              // Else, open file and append output.
        file_1.open( "time.txt" );
    }
    else
    {
        file_1.open( "time.txt", std::ios::app );
    }
    
    if ( !check_file_2.is_open() )
    {
        file_2.open( "virus.txt" );
    }
    else
    {
        file_2.open( "virus.txt", std::ios::app );
    }
    
    if ( !check_file_3.is_open() )
    {
        file_3.open( "gen.txt" );
    }
    else
    {
        file_3.open( "gen.txt", std::ios::app );
    }
    
    if ( !check_file_4.is_open() )
    {
        file_4.open( "parent.txt" );
    }
    else
    {
        file_4.open( "parent.txt", std::ios::app );
    }
    
    // Loops through vector database of individuals and outputs appropriate
    // data to the output files.
    for (std::vector<individual*>::iterator it = population.begin(); it != population.end(); ++it)
    {
        
        file_1 << std::setprecision(4) << (*it)->t_0 << '\t';

        for (int j = 0; j < K-1; j++)
        {
            file_2 << (*it)->v_0[j] << ", ";
            
        }
        
        file_2 << (*it)->v_0[K-1] << std::endl;
        file_3 << (*it)->gen << '\t';

        
        if ((*it)->parent)
            file_4 << (*it)->parent->key << '\t';
    }
    
    // Since we're done outputting data, close the files and return to calling
    // function.
    file_1.close();
    file_2.close();
    file_3.close();
    file_4.close();
}

// Simply formats and outputs data to file.
void CrossScaleSim::output_formatted()
{
    std::ifstream check_file_formatted ("formatted.txt");
    std::ofstream formatted_output;
    
    if ( !check_file_formatted.is_open() ) // If the file doesn't exist, create it.
    {                              // Else, open file and append output.
        formatted_output.open( "formatted.txt" );
    }
    else
    {
        formatted_output.open( "formatted.txt", std::ios::app );
    }
    
    for (std::vector<individual*>::iterator it = population.begin() ; it != population.end(); ++it)
    {
        
        
        for (int j = 0; j < K; j++)
        {
            formatted_output << (*it)->v_0[j] << ", ";
            
        }

        formatted_output << (*it)->t_0 << ", ";
        formatted_output << (*it)->key << ", ";
        
        if ((*it)->parent)
            formatted_output << (*it)->parent->key << ", ";
        else
            formatted_output << "00, ";
        
        formatted_output << (*it)->gen << std::endl;
        
    }
    
    formatted_output.close();
}

// Simply formats and outputs data to console.
void CrossScaleSim::output_map()
{
    std::cout << '|' << std::setw(2) << " Viral Type: ";
    std::cout << std::setfill(' ') << "[1] " << '\t' << "[2] " << '\t' ;
    std::cout << '|' << std::setw(3) << " Time " << std::setw(5) << '|' << " Key ";
    std::cout << '|' << std::setw(5) << " Parent " << '|' << " Gen " << '|' << std::endl;
    
    for (std::vector<individual*>::iterator it = population.begin() ; it != population.end(); ++it)
    {
        
        std::cout << std::setw(8) << '\t';
        
        for (int j = 0; j < K; j++)
        {
            std::cout << std::setw(8) << std::setfill(' ') << (*it)->v_0[j];
            
        }
        std::cout << '\t' << std::setw(8) << std::setprecision(4) << std::setfill(' ') << (*it)->t_0;
        std::cout << std::setw(6) << std::setfill(' ') << (*it)->key;
        
        if ((*it)->parent)
            std::cout << std::setw(6) << std::setfill(' ') << '\t' << (*it)->parent->key;
        else
            std::cout << std::setw(6) << std::setfill(' ') << '\t' << "NULL";
        
        std::cout << std::setw(4) << std::setfill(' ') << '\t' << (*it)->gen;
        std::cout << std::endl;
        
    }
}

std::vector<double> CrossScaleSim::betweenHost( double time_i, double time_f, std::vector<int>types )
{
    ODESolver VPop(types, time_i, time_f);
    std::vector<double> wtypes = VPop.Iterate();
    
    return wtypes;
}

std::vector<int> CrossScaleSim::withinHost( individual* parent )
{
    std::vector<int> types = parent->v_0;
    ODESolver initial(types, T_0, T);
    types = initial.Iterate( 1.0 );
    parent->t_crit = parent->t_0+initial.get_time();
    parent->t_abs = initial.get_time();
    
    return types;
}

void CrossScaleSim::transmission( std::vector<int> parentvir, individual* parent )
{
    double t_crit = parent->t_crit;
    int parent_gen = parent->gen;
    int contacts = poiss(gen); // Generates number from poisson distro.
    std::vector<double> times;
    
    //Uncomment the line below to debug.
    //std::cout << "Total # of contacts: " << contacts << std::endl;
    
    for (int i = 0; i < contacts; i++)
    {
        times.push_back(unif(gen)); // Generates number from uniform distro.
    }

    for (int i = 0; i < times.size(); i++)
    {
        std::vector<double> w_temp = betweenHost(0.0, times[i], parentvir);
        
        std::vector<int> wtypes;
        
        for (int k = 0; k < K; k++)
        {
            //unsigned long index = rand() % 5;
            //lambda = relative transmissibility*mean viral load
            double lambda = B[k]*w_temp[k];
            
            //Uncomment the line below to debug.
            //std::cout << "Lamdba is equal to:" << lambda << std::endl;
            
            std::poisson_distribution<> p_temp(lambda); // Create temp poisson distro.
            int particles = p_temp(gen); // Generate poisson number.
            
            //Uncomment the line below to debug.
            //std::cout << "# of particles from virus type " << k << " is: " << particles << std::endl;
            
            wtypes.push_back(particles);
        }
        
        
        int total_particles = 0;
        
        for (std::vector<int>::iterator it = wtypes.begin() ; it != wtypes.end(); ++it)
        {
            total_particles += *it;
        }

        if ( total_particles == 0 )
        {
            return;
        }
        
        else
        { // Push new individuals into database.
            (*parent).children.push_back( new individual( ( t_crit+times[i]), wtypes, parent, (population.size()+children.size()+1), parent_gen+1 ) );
            children.push(parent->children.back());
        }
    }
}

void CrossScaleSim::avgTimeToCrit()
{    
    for (std::vector<individual*>::iterator it = population.begin() ; it != population.end(); ++it)
    {
        t_crit.push_back((*it)->t_abs);
    }
}
