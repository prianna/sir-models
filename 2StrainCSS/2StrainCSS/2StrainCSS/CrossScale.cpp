//
//  CrossScale.cpp
//  Source code for CrossScaleSim class.

#include "CrossScale.h"


CrossScaleSim::CrossScaleSim()
: generations(0)
{
    std::random_device rd;
    randomgen.seed(rd()); // Seed Mersenne Twister pseduo RNG with random device.
	V_0 = VPOP;
    population.push_back(new individual(0.0, V_0, NULL, 1, generations));
    individual *parent = population.front();
    ODESolver initial(parent->v_0, T_0, T_max);
    std::vector<int> v_crit = initial.Iterate( T_max );
    parent->t_crit = initial.get_time();
    parent->t_abs = initial.get_time();
    transmission(v_crit, parent);
}

// Main iterative loop. This loop cycles through the a queue of infected
// individuals, checking whether they have reached the critical threshold for
// infection (V_CRIT). These individuals then enter the transmission phase
// or are allowed to continue incubating the virus until V_CRIT is reached.
// Calls withinHost() to carry out within-host incubation period and
// the transmission() function to calculate bottleneck width.
void CrossScaleSim::cycle()
{
    // Checking to make sure that our vector of individuals
    // doesn't grow too large.
    if (population.size() >= POP_MAX)
    {
        return;
    }
    
    // If there are newly infected individuals in the queue, waiting to transmit,
    // we check their viral load to make sure they have reached V_CRIT.
    else if ( !children.empty() )
    {
        while ( !children.empty() && generations <= MAX_GENS && population.size() < POP_MAX )
        {
            individual *parent = children.front();
            generations = parent->gen;
            
            if( generations < MAX_GENS )
            {
                double v_load = viral_Load(parent->v_0);
                
                if ( v_load < V_CRIT )
                {
                    std::vector<int> types = withinHost(parent);
                    transmission(types, parent);
                }
                else
                {
                    transmission(parent->v_0, parent);
                }
            }
            
            // Adds current child to vector of parents.
            population.push_back(parent);
            
            // Dequeues current child from the list of children waiting to
            // transmit.
            children.pop();
        }
        return;
    }
}


// Calls ODE integrator to run from T until T_max or until V_crit is attained.
// Called by cycle()
std::vector<int> CrossScaleSim::withinHost( individual* parent )
{
    std::vector<int> types = parent->v_0;
    ODESolver initial(types, T_0, T_max);
    
    types = initial.Iterate( T_max );
    
    parent->t_abs = initial.get_time();
    parent->t_crit = (parent->t_0)+(parent->t_abs);
    
    return types;
}


// Calls ODE integrator to run from T until time specified by transmission function.
// Called by transmission()
std::vector<double> CrossScaleSim::betweenHost( double time_i, double time_f, std::vector<int> types )
{
    ODESolver VPop(types, time_i, time_f);
    std::vector<double> wtypes = VPop.Iterate();
    
    // Uncomment below to debug.
    /*for (std::vector<double>::iterator it = wtypes.begin() ; it != wtypes.end(); ++it)
    {
        std::cout << "# of particles is: " << *it << std::endl;

    }*/
    
    return wtypes;
}


// Decides when transmission happens, how many particles are transmitted, and
// how many people are infected.
// Calls betweenHost() to simulate growth.
// Modifies parent's vector of children created and enqueues new children, who will
// enter their own transmission cycle.
void CrossScaleSim::transmission( std::vector<int> parentvir, individual* parent )
{
    // Poisson distribution for # contacts.
    std::poisson_distribution<int> poiss(T*BETA_0);
    int contacts = poiss(randomgen);
    parent->contacts = contacts;
    
    //Uncomment the line below to debug.
    //std::cout << "Total # of contacts: " << contacts << std::endl;
    
    // Uniform distribution for transmission times.
    std::uniform_real_distribution<double> unif(0.0,T);
    std::vector<double> times;

    for (int i = 0; i < contacts; i++)
    {
        // Generates time from uniform distro for each potential contact.
        times.push_back(unif(randomgen));
    }
    
    //Uncomment the line below to debug.
    //std::cout << "Total number of times: " << times.size() << std::endl;

    //Iterating through vector of times.
    for (std::vector<double>::iterator it = times.begin() ; it != times.end(); ++it)
    {
        // Creates a temporary vector to hold the raw viral load for the potential child.
        std::vector<double> w_temp = betweenHost(0.0, (*it), parentvir);
        
        // This will be the potential new child with poisson-normalized load.
        std::vector<int> wtypes;
        
        int total_particles = 0;

        for (int k = 0; k < K; k++)
        {
            // lambda := product of relative transmissibility and mean viral load.
            double lambda = B[k]*w_temp[k];
            
            // Create temp poisson distro for determination of bottleneck width.
            std::poisson_distribution<int> p_temp(lambda);
            int particles = p_temp(randomgen);
            // Keeping track of the total # of particles transmitted.
            total_particles += particles;
            
            // Populate child's viral load vector with particles.
            wtypes.push_back(particles);

            //Uncomment the line below to debug.
            //std::cout << "# of particles from virus type " << k << " is: " << particles << std::endl;
        }
        
        // If there was successful transmission of at least 1 particle.
        if ( total_particles > 0 )
        {
            // Push new individuals into parent's personal children database.
            parent->children.push_back( new individual( (*it), wtypes, parent, (population.size()+children.size()+1), (parent->gen)+1 ) );
            // Push new individuals into the queue of children waiting to transmit.
            children.push(parent->children.back());
        }
    }
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
        
        formatted_output << (*it)->gen << ", " << (*it)->children.size() <<
        ", " << (*it)->contacts <<std::endl;
        
    }
    
    formatted_output.close();
}

// Simply formats and outputs data to console.
void CrossScaleSim::output_map()
{
    std::cout << '|' << std::setw(2) << " Viral Type: ";
    std::cout << std::setfill(' ') << "[1] " << '\t' << "[2] " << '\t' ;
    std::cout << '|' << std::setw(3) << " Time " << std::setw(5) << '|' << " Key ";
    std::cout << '|' << std::setw(5) << " Parent " << '|' << " Gen " << '|'
    << " Children " << '|' <<std::endl;
    
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
        std::cout << std::setw(4) << std::setfill(' ') << '\t' << (*it)->children.size();
        
        std::cout << std::endl;
        
    }
}


void CrossScaleSim::avgTimeToCrit()
{    
    for (std::vector<individual*>::iterator it = population.begin() ; it != population.end(); ++it)
    {
        t_crit.push_back((*it)->t_abs);
    }
}
