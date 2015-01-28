/* ****************************************************************

This is the C version of program 2.2 from page 27 of
"Modeling Infectious Disease in humans and animals"
by Keeling & Rohani.

It is the simple SIR epidemic with equal births and deaths.

This code is written to be simple, transparent and readily compiled.
Far more elegant and efficient code can be written.

This code can be compiled using gnu or intel C compilers:
icc -o Program_2_2  Program_2_2.c  -lm
g++ -o Program_2_2  Program_2_2.c  -lm

If you have any difficulties, we suggest you set DEBUG to 1 or -1

******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctime>


#define DEBUG 1

const int NUM_EVENTS = 5;

// Set up basic parameters
double beta=1; // Transmission term
double gamm=1.0/10; // Mutation rate for WT virions
double mu=1.0; // Death rate & birth rate
int S0= 40; // Uninfected cells
int I0= 5; // WT infected cells
int IM0= 0; // Mutated infected cells
double MaxTime=5.0;

// Set up variables and rates of change
double t;
int S,I,IM;


void Read_in_parameters(FILE *);
void Perform_Checks();
void Output_Data(FILE *);
double Gillespie();


// The main program

int main(int argc, char** argv)
{

  srand(time(0));

  char FileName[1000],str[200];
  double Every;
  FILE *in, *out;
   
  /* Find the given parameter file name and open the file 
     otherwise use parameters set at top of program */
  if(argc>1)
    {
      strcpy(FileName, *(argv+1));  
      
      if((in = fopen(FileName,"r"))==NULL)
	{
	  printf("Cannot read file %s\n",FileName);
	  exit(1);
	}

      /* Read in the parameters */
      
      Read_in_parameters(in);
    
      if(DEBUG) printf("\nReading parameters from file %s\n",FileName);
    }
  else
    {
      if(DEBUG) printf("\nUsing default parameters\n");
    }

  if(DEBUG) printf("beta=%g\ngamma=%g\nmu=%g\nInitial S=%d\nInitial I=%d\nRuns until time %g\n\n",beta,gamm,mu,S0,I0,MaxTime);

  /* Check all parameters are OK & set up intitial conditions */

  Perform_Checks();
  
  S=S0; I=I0; IM=gamm*I0;


  if((out = fopen("Output","w"))==NULL)
    {
      printf("Cannot open file Output to write data\n");
      exit(1);
    }

  /* The main iteration routine */

  t=0;

  Output_Data(out);

  do
    {
	  double dt = Gillespie();
      t+=dt;
      /* If time has moved on sufficiently, output the current data */
      //if( floor(t/Every) > floor((t-dt)/Every) )
	//{
	  Output_Data(out);
	//}
      //else
	//{
	  if(DEBUG<0) Output_Data(stdout);
	//}
      
    }
  while(t<5.0);

  Output_Data(out);

  fclose(out);
}



void Read_in_parameters(FILE *in)
{
  char str[200];
  
  fscanf(in,"%s",str);
  fscanf(in,"%lf",&beta);
  
  fscanf(in,"%s",str);
  fscanf(in,"%lf",&gamm);
  
  fscanf(in,"%s",str);
  fscanf(in,"%lf",&mu);
  
  fscanf(in,"%s",str);
  fscanf(in,"%lf",&S0);
  
  fscanf(in,"%s",str);
  fscanf(in,"%lf",&I0);
  
  fscanf(in,"%s",str);
  fscanf(in,"%lf",&MaxTime);
  
  fclose(in);
}



void Perform_Checks()
{
  if(S0<=0) {printf("ERROR: Initial level of susceptibles (%g) is less than or equal to zero\n",S0); exit(1);}

  if(I0<=0) {printf("ERROR: Initial level of infecteds (%g) is less than or equal to zero\n",I0); exit(1);}

  if(beta<=0) {printf("ERROR: Transmission rate beta (%g) is less than or equal to zero\n",beta); exit(1);}

 if(gamm<=0) {printf("ERROR: Recovery rate gamma (%g) is less than or equal to zero\n",gamm); exit(1);}

 if(mu<0) {printf("ERROR: Birth / Death rate mu (%g) is less than zero\n",mu); exit(1);}

  if(MaxTime<=0) {printf("ERROR: Maximum run time (%g) is less than or equal to zero\n",MaxTime); exit(1);}

  if(beta<gamm+mu) {printf("WARNING: Basic reproductive ratio (R_0=%g) is less than one\n",beta/(gamm+mu));}
}


void Output_Data(FILE *out)
{
  
  if(out!=stdout) fprintf(out,"%g \t %d \t %d \t %d \n",t,S,I,IM);

  if(DEBUG) printf("%g \t %d \t %d \t %d \n",t,S,I,IM);
}

double Gillespie()
{
	double rates[NUM_EVENTS] = { beta*S*I, beta*S*IM, gamm*beta*S*I, mu*I, mu*IM };
	//double dt[NUM_EVENTS];
	double totalRate = rates[1]+rates[2]+rates[3]+rates[4]+rates[5];
	double r = (((double) rand() / (RAND_MAX)) + 1);
	double dt = (-1/totalRate)*log( r );
	//double last = 0;
	double p = (((double) rand() / (RAND_MAX)) + 1);



	for( int i = 0; i < NUM_EVENTS; i++ )
	{
		if( p < rates[i]/totalRate )
		{
			switch (i)
			{
			case 0:
				S--;
				I++;
				break;

			case 1:
				S--;
				IM++;
				break;

			case 2:
				IM++;
				I--;
				break;

			case 3:
				I--;
				break;

			case 4:
				IM--;
				break;

			default:
				break;
			}
		}

		//last += (rates[i]);
	}

	return dt;

	//for( int i = 0; i < NUM_EVENTS; i++ )
	//{
		//dt[i] = (log( (float)rand()/RAND_MAX ))*(-1/totalRate);
	//}

}
