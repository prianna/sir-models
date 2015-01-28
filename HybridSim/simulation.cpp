//
//  simulations.cpp
//  HybridSim
//
//  Created by Prianna Ahsan on 5/27/13.
//  Copyright (c) 2013 Prianna Ahsan. All rights reserved.
//

#include "simulation.h"

void Simulation::Output_Data(FILE *output)
{
    if(output!=stdout)
        fprintf(output,"%g \t %g \t %g \t %g \n",t,S,I1,I2);
    
    if(DEBUG)
        printf("%g \t %g \t %g \t %g \n",t,S,I1,I2);
}

void Simulation::Iterate( FILE *output )
{
    Output_Data(output);
    while (t < MAX_TIME )
    {
        if (I1 < 1000 || I2 < 1000)
        {
            dt = Gillespie();
            Runge_Kutta(dt);
            t+=dt;
            Output_Data(output);
        }
        else
        {
            int n = (I1 >=I2 ? I1 : I2);
            Runge_Kutta(dt, n);
            t+=dt;
            Output_Data(output);
        }
    }
    Output_Data(output);
    return;
}

void Simulation::Diff(double Pop[3])
{
// Set up temporary variables to make the equations look neater
        
    double tmpS, tmpI1, tmpI2;
    
    tmpS=Pop[0]; tmpI1=Pop[1]; tmpI2=Pop[2];
    
    /* The differential equations */
    
    dPop[0] = lambda-dr*tmpS-beta1*tmpS*tmpI1-beta2*tmpS*tmpI2; // dS/dt
    dPop[1] = (1-mu)*beta1*tmpS*tmpI1-del*tmpI1;  // dI1/dt
    dPop[2] = mu*beta1*tmpS*tmpI1+beta2*tmpS*tmpI2-del*tmpI2; // dI2/dt
    
    return;
}

void Simulation::Runge_Kutta(double step)
{
    double dPop1[3], dPop2[3], dPop3[3], dPop4[3];
    double tmpPop[3], initialPop[3];
    
    /* Step sizes causing issues */
    
    initialPop[0]=S; initialPop[1]=I1; initialPop[2]=I2;
    tmpPop[0]=initialPop[0];
    tmpPop[1]=initialPop[1];
    tmpPop[2]=initialPop[2];
    
    Diff(initialPop);

    dPop1[0]=dPop[0];
    tmpPop[0]=initialPop[0]+step*dPop1[0]/2;
    
    Diff(tmpPop);
    dPop2[0]=dPop[0];
    tmpPop[0]=initialPop[0]+step*dPop2[0]/2;
    
    Diff(tmpPop);
    dPop3[0]=dPop[0];
    tmpPop[0]=initialPop[0]+step*dPop3[0];
    
    Diff(tmpPop);
    dPop4[0]=dPop[0];
    tmpPop[0]=initialPop[0]+(dPop1[0]/6 + dPop2[0]/3 + dPop3[0]/3 + dPop4[0]/6)*step;
    
    
    S=tmpPop[0];
    
    return;
}

void Simulation::Runge_Kutta(double step, int n)
{
    if (n >= 1000 )
    {
        step=0.1/((beta1+beta2+del+mu)*S0);
        int i;
        double dPop1[3], dPop2[3], dPop3[3], dPop4[3];
        double tmpPop[3], initialPop[3];
        
        /* Integrates the equations one step, using Runge-Kutta 4
         Note: we work with arrays rather than variables to make the
         coding easier */
        
        initialPop[0]=S; initialPop[1]=I1; initialPop[2]=I2;
        
        Diff(initialPop);
        for(i=0;i<3;i++)
        {
            dPop1[i]=dPop[i];
            tmpPop[i]=initialPop[i]+step*dPop1[i]/2;
        }
        
        Diff(tmpPop);
        for(i=0;i<3;i++)
        {
            dPop2[i]=dPop[i];
            tmpPop[i]=initialPop[i]+step*dPop2[i]/2;
        }
        
        Diff(tmpPop);
        for(i=0;i<3;i++)
        {
            dPop3[i]=dPop[i];
            tmpPop[i]=initialPop[i]+step*dPop3[i];
        }
        
        Diff(tmpPop);
        for(i=0;i<3;i++)
        {
            dPop4[i]=dPop[i];
            
            tmpPop[i]=initialPop[i]+(dPop1[i]/6 + dPop2[i]/3 + dPop3[i]/3 + dPop4[i]/6)*step;
        }
        
        
        S=tmpPop[0]; I1=tmpPop[1]; I2=tmpPop[2];
        
        return;
    }
    else
    {
        Runge_Kutta(step);
        return;
    }
}

double Simulation::Gillespie()
{
	double rates[NUM_EVENTS] = {(1-mu)*beta1*S*I1, del*I1, mu*beta1*S*I1, beta2*S*I2, del*I2};
	double totalRate = (rates[0]+rates[1]+rates[2]+rates[3]+rates[4]);
	double r = (((double) rand() / (RAND_MAX)));
	double dt = (-1/totalRate)*log( r );
	double last = 0;
    double p = (((double) rand() / (RAND_MAX)));
    
    
	for( int i = 0; i < NUM_EVENTS; i++ )
	{
		if( (last/totalRate) < p && p < (last+rates[i])/totalRate )
            
		{
			switch (i)
			{
                case 0:
                    I1++;
                    return dt;
                    
                case 1:
                    I1--;
                    return dt;
                    
                case 2:
                    I2++;
                    return dt;
                    
                case 3:
                    I2++;
                    return dt;
                    
                case 4:
                    I2--;
                    return dt;
                    
                default:
                    break;
			}
		}
        
		last += rates[i];
	}
    
	return dt;
}
