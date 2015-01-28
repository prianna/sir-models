//
//  main.cpp
//  HybridSim
//
//  Created by Prianna Ahsan on 5/27/13.
//  Copyright (c) 2013 Prianna Ahsan. All rights reserved.
//


#include "simulation.h"
#include "globals.h"

#define DEBUG 1

void Read_in_parameters(FILE *in);

int main(int argc, const char * argv[])
{
    srand(time(0));
    
    for (int i = 0; i < MAX_RUNS; i++)
    {
        char FileName[1000];
        FILE *input, *output;
        std::string out;
        Simulation sim;
        
        /* Find the given parameter file name and open the file
         otherwise use parameters set at top of program */
        if(argc>1)
        {
            strcpy(FileName, *(argv+1));
            
            if((input = fopen(FileName,"r"))==NULL)
            {
                printf("Cannot read file %s\n",FileName);
                exit(1);
            }
            
            Read_in_parameters(input);
            
            if(DEBUG)
                printf("\nReading parameters from file %s\n",FileName);
        }
        
        else
        {
            if(DEBUG)
                printf("\nUsing default parameters\n");
        }
        
        if(DEBUG)
        {
            printf("beta1=%g\nbeta2=%g\nmu=%g\nlambda=%g\ndr=%g\ndel=%g\nInitial S=%g\nInitial I1=%g\nInitial I2=%g\nRuns until time %g\n\n",beta1,beta2,mu,lambda, dr, del, S0,I10, I20, MAX_TIME);

            std::stringstream ss;
            ss << "Output" << i;
            out = ss.str();
        }
        
        
        if((output = fopen(out.c_str(),"w"))==NULL)
        {
            printf("Cannot open file to write data\n");
            exit(1);
        }
        
        sim.Iterate(output);
        fclose(output);
    }
    
}




void Read_in_parameters(FILE *in)
{
    char str[200];
    
    fscanf(in,"%s",str);
    fscanf(in,"%lf",&beta1);
    
    fscanf(in,"%s",str);
    fscanf(in,"%lf",&beta2);
    
    fscanf(in,"%s",str);
    fscanf(in,"%lf",&mu);
    
    fscanf(in,"%s",str);
    fscanf(in,"%lf",&lambda);
    
    fscanf(in,"%s",str);
    fscanf(in,"%lf",&dr);
    
    fscanf(in,"%s",str);
    fscanf(in,"%lf",&del);
    
    fscanf(in,"%s",str);
    fscanf(in,"%lf",&S0);
    
    fscanf(in,"%s",str);
    fscanf(in,"%lf",&I10);
    
    fscanf(in,"%s",str);
    fscanf(in,"%lf",&I20);
    
    fscanf(in,"%s",str);
    fscanf(in,"%lf",&MAX_TIME);
    
    fclose(in);
}

