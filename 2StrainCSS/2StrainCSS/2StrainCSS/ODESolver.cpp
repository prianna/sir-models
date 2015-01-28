#include "ODESolver.h"
#include "globals.h"

using namespace boost::numeric::odeint;

typedef controlled_runge_kutta< runge_kutta_fehlberg78< state_type > > stepper_type;

void System( state_type &V , state_type &dVdt , double t )
{
    dVdt[0] = RATES[0]*V[0] + (V[1]-V[0])*mu; // dWT/dt
    //dVdt[1] = RATES[1]*V[1] + (V[0]-V[1])*mu;  // dM/dt
}

ODESolver::ODESolver(std::vector<int> init_types, double t_init, double t_final)
:t(t_init), t_final(t_final), step(0.01)
{
    for (int i = 0; i < K; i++)
    {
        V.push_back(init_types[i]);
    }
}

ODESolver::ODESolver()
:t(T_0), t_final(T), step(0.01)
{
}

std::vector<double> ODESolver::Iterate()
{
    
    prepare_output( formatted_output, "BetweenHost.txt" );
    integrate_adaptive( stepper_type() , System , V , t , t_final , step,
                                streaming_observer( formatted_output ));
    
    end_output( formatted_output );
    std::vector<double> wtypes;
    
    double v_total = 0.0;
    
    for (int i = 0; i < K; i++)
    {
        v_total += V[i];
    }
    
    if (v_total < N)
    {
        for (int i = 0; i < K; i++)
        {
            wtypes.push_back(V[i]);
        }
        return wtypes;
    }
    else
    {
        for (int i = 0; i < K; i++)
        {
            double tempV = ( (V[i] / v_total) * N );
            wtypes.push_back(tempV);
        }
        return wtypes;
    }
}

std::vector<int> ODESolver::Iterate( double T )
{
    //runge_kutta_fehlberg78< state_type > rkf_78;
    stepper_type stepper;
    std::ofstream output;
    prepare_output( output, "WithinHost.txt" );
    streaming_observer out(output);
    double dt = 0.01;
    
    double v_total = 0.0;
    
    for (int i = 0; i < K; i++)
    {
        v_total += V[i];
    }
    
    while ( v_total < V_CRIT )
    {
        stepper.try_step( System , V , t_final , dt );
        
        v_total = 0.0;
        for (int i = 0; i < K; i++)
        {
            v_total += V[i];
        }
        
        out(V, t);
    }
    
    /*do
    {
        stepper.try_step( System , V , t_final , dt );
        t+= dt;
        
        
        v_total = 0.0;
        for (int i = 0; i < K; i++)
        {
            v_total += V[i];
        }
        
        out(V, t);
    }
    while ( v_total < V_CRIT );*/
    
    std::vector<int> wtypes;

    for (int i = 0; i < K; i++)
    {
        int W = round(V[i]);
        wtypes.push_back(W);
    }
    
    return wtypes;
}

std::vector<double> ODESolver::IterateR( double T0, double TF ) const
{
    std::vector<double> types = V;
    
    integrate_adaptive( stepper_type() , System , types , T0 , TF , step);
    
    double v_total = 0.0;
    
    for (int i = 0; i < K; i++)
    {
        v_total += types[i];
    }
    
    if (v_total < N)
        return types;
    
    else
    {
        for (int i = 0; i < K; i++)
        {
            double tempV = ( (types[i] / v_total) * N );
            types[i] = tempV;
        }
        return types;
    }
}

void ODESolver::prepare_output( std::ofstream &output, std::string filename )
{
    std::ifstream check_file_formatted (filename);
    
    if ( !check_file_formatted.is_open() ) // If the file doesn't exist, create it.
    {                              // Else, open file and append output.
        output.open( filename );
    }
    else
    {
        output.open( filename, std::ios::app );
    }
}

/*
void ODESolver::Diff(double Pop[K])
{
    // Set up temporary variables to make the equations look neater
    
    double tmpV1, tmpV2;
    
    tmpV1=Pop[0]; tmpV2=Pop[1];
    
    
    dPop[0] = RATES[0]*V[0] + (V[1]-V[0])*mu; // dWT/dt
    dPop[1] = RATES[1]*V[1] + (V[0]-V[1])*mu;  // dM/dt
    
    return;
}


void ODESolver::Runge_Kutta(double step)
{
    int i;
    double dPop1[K], dPop2[K], dPop3[K], dPop4[K];
    double tmpPop[K], initialPop[K];
    

    
    initialPop[0]=V[0]; initialPop[1]=V[1];
    
    //K1
    Diff(initialPop);
    
    //K2
    for(i=0;i<K;i++)
    {
        dPop1[i]=dPop[i];
        tmpPop[i]=initialPop[i]+step*dPop1[i]/2;
    }
    Diff(tmpPop);
    
    //K3
    for(i=0;i<K;i++)
    {
        dPop2[i]=dPop[i];
        tmpPop[i]=initialPop[i]+step*dPop2[i]/2;
    }
    Diff(tmpPop);
    
    //K4
    for(i=0;i<K;i++)
    {
        dPop3[i]=dPop[i];
        tmpPop[i]=initialPop[i]+step*dPop3[i];
    }
    Diff(tmpPop);
    
    //Y_(i+1)
    for(i=0;i<K;i++)
    {
        dPop4[i]=dPop[i];
        tmpPop[i]=initialPop[i]+(dPop1[i]/6 + dPop2[i]/3 + dPop3[i]/3 + dPop4[i]/6)*step;
    }
    
    
    V.at(0)=tmpPop[0]; V.at(1)=tmpPop[1];
    
    return;
}


void ODESolver::Runge_Kutta_Fehlberg(double h)
{
    int i;
    double dPopT[6][3];
    double tmpPop[K], initialPop[K];
    
    
    initialPop[0]=V[0]; initialPop[1]=V[1];
    
    //K1
    Diff(initialPop);
    
    //K2
    for(i=0;i<K;i++)
    {
        dPopT[0][i]=dPop[i];
        tmpPop[i]=initialPop[i]+step*dPopT[0][i]/2;
    }
    Diff(tmpPop);
    
    //K3
    for(i=0;i<K;i++)
    {
        dPopT[1][i]=dPop[i];
        tmpPop[i]=initialPop[i]+step*dPopT[1][i]/2;
    }
    Diff(tmpPop);
    
    //K4
    for(i=0;i<K;i++)
    {
        dPopT[2][i]=dPop[i];
        tmpPop[i]=initialPop[i]+step*dPopT[2][i];
    }
    Diff(tmpPop);
    
    //Y_(i+1)
    for(i=0;i<K;i++)
    {
        dPopT[3][i]=dPop[i];
        tmpPop[i]=initialPop[i]+(dPopT[0][i]/6 + dPopT[1][i]/3 + dPopT[2][i]/3 + dPopT[3][i]/6)*step;
    }
    
    
    V.at(0)=tmpPop[0]; V.at(1)=tmpPop[1];
    
    return;
}
*/

