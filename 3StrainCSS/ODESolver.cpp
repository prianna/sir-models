#include "ODESolver.h"
#include "globals.h"

ODESolver::ODESolver(std::vector<int> init_types, double t_init, double t_final)
:t(t_init), t_final(t_final)
{
    for (int i = 0; i < K; i++)
    {
        V.push_back(init_types[i]);
        dPop[i]=0;
    }
    
    step=0.00001/((RATES[0]+RATES[1]+RATES[2])*mu);
}

std::vector<double> ODESolver::Iterate()
{
    do
    {
        Runge_Kutta(step);
        t+=step;
    }
    while(t<t_final);
    
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
            double tempV = ((V[i]/v_total)*N);
            wtypes.push_back(tempV);
        }
        return wtypes;
    }
}

std::vector<int> ODESolver::Iterate( double total )
{
    double v_total;
    
    do
    {
        Runge_Kutta(step);
        t+=step;
        v_total = 0.0;
        for (int i = 0; i < K; i++)
        {
            v_total += V[i];
            
            //if (total == 2.0)
              //  std::cout << "Type " << i+1 << ": " << V[i] << " ";
        }
        
       // if ( total == 2.0 )
            //std::cout << std::endl << "Time: " << t << std::endl;
    }
    while(v_total < V_CRIT);
    
    std::vector<int> wtypes;

    for (int i = 0; i < K; i++)
    {
        int W = round(V[i]);
        wtypes.push_back(W);
    }
    
    return wtypes;
}


void ODESolver::Diff(double Pop[K])
{
    // Set up temporary variables to make the equations look neater
    
    double tmpV1, tmpV2, tmpV3;
    
    tmpV1=Pop[0]; tmpV2=Pop[1]; tmpV3=Pop[2];
    
    /* The differential equations */
    
    dPop[0] = RATES[0]*V[0] + mu*(V[1]-V[0]) + (mu*mu)*(V[2]-V[0]); // dV1/dt
    dPop[1] = RATES[1]*V[1] + mu*(V[0] + V[2] - 2*V[1]);  // dV2/dt
    dPop[2] = RATES[2]*V[2] + mu*(V[1] - V[2]) + (mu*mu)*(V[0]-V[2]); //dV3/dt
    
    return;
}

void ODESolver::Runge_Kutta(double step)
{
    int i;
    double dPop1[K], dPop2[K], dPop3[K], dPop4[K];
    double tmpPop[K], initialPop[K];
    
    /* Integrates the equations one step, using Runge-Kutta 4
     Note: we work with arrays rather than variables to make the
     coding easier */
    
    initialPop[0]=V[0]; initialPop[1]=V[1]; initialPop[2]=V[2];
    
    Diff(initialPop);
    for(i=0;i<K;i++)
    {
        dPop1[i]=dPop[i];
        tmpPop[i]=initialPop[i]+step*dPop1[i]/2;
    }
    
    Diff(tmpPop);
    for(i=0;i<K;i++)
    {
        dPop2[i]=dPop[i];
        tmpPop[i]=initialPop[i]+step*dPop2[i]/2;
    }
    
    Diff(tmpPop);
    for(i=0;i<K;i++)
    {
        dPop3[i]=dPop[i];
        tmpPop[i]=initialPop[i]+step*dPop3[i];
    }
    
    Diff(tmpPop);
    for(i=0;i<K;i++)
    {
        dPop4[i]=dPop[i];
        
        tmpPop[i]=initialPop[i]+(dPop1[i]/6 + dPop2[i]/3 + dPop3[i]/3 + dPop4[i]/6)*step;
    }
    
    
    V.at(0)=tmpPop[0]; V.at(1)=tmpPop[1]; V.at(2)=tmpPop[2];
    
    return;
}


