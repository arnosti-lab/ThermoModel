/* --------------------------------------------------------- */
/* --------------- A Very Short Example -------------------- */
/* --------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h> /* free() */
#include <math.h>
#include "cmaes_interface.h"
#include "ThermoModel.h"
#include "All_Inputs.h"

double fitfun(double const *x, int dim); 

/* the objective (fitness) function to be minized */
double fitfun(double const *x, int N) { /* function "cigtab" */
	double *pos = new double[N];
	int i;
	int typeC = 2, numC = 12;
	int typeQ = 1, numQ = 8;

    double lowerSF = 0;
	double lowerC = 0;
	
	for(i=0; i<N; i++)
		pos[i] = x[i]*x[i];
	//puts lower bound on SCALING FACTORS
	for(i=0; i<3; i++)
		pos[i] = pos[i]+lowerSF;
	//puts lower bound on COOPERATIVITIES
	for(i=3; i<3+numC; i++)
		pos[i] = pos[i]+lowerC;
	//if binned quenching....it will set your quenching parameters in [0,1]
	if(typeQ==1) {
		for (i=N-numQ; i<N; i++)
			pos[i] = exp(-pos[i]);
	}
	double SSE = runobjective(pos,N,typeC,typeQ);
	delete[] pos;
	return SSE;
}

/* the optimization loop */
int main(int argn, char **args) {
	initialize_thermo();
    cmaes_t evo; /* an CMA-ES type struct or "object" */
    double *arFunvals, *const*pop, *xfinal;

    /* Initialize everything into the struct evo, 0 means default */
    arFunvals = cmaes_init(&evo, 0, NULL, NULL, 0, 0, "initials.par"); 
    printf("%s\n", cmaes_SayHello(&evo));
    cmaes_ReadSignals(&evo, "signals.par");  /* write header and initial values */

    /* Iterate until stop criterion holds */
    while(!cmaes_TestForTermination(&evo)) { 
        /* generate lambda new search points, sample population */
        pop = cmaes_SamplePopulation(&evo); /* do not change content of pop */

        /* Here you may resample each solution point pop[i] until it
	    becomes feasible, e.g. for box constraints (variable
	    boundaries). function is_feasible(...) needs to be
	    user-defined.  
	    Assumptions: the feasible domain is convex, the optimum is
	    not on (or very close to) the domain boundary, initialX is
	    feasible and initialStandardDeviations are sufficiently small
	    to prevent quasi-infinite looping.
        */
        /* for (int i = 0; i < cmaes_Get(&evo, "popsize"); ++i) while (!is_feasible(pop[i])) cmaes_ReSampleSingle(&evo, i); 
        */

        int lambda = cmaes_Get(&evo, "lambda");
        int dimn = (int) cmaes_Get(&evo, "dim");

        /* evaluate the new search points using fitfun from above */
        for (int i = 0; i < lambda; i++) {
	        arFunvals[i] = fitfun(pop[i], dimn);
        }

        /* update the search distribution used for cmaes_SampleDistribution() */
        cmaes_UpdateDistribution(&evo, arFunvals);  

        /* read instructions for printing output or changing termination conditions */ 
        cmaes_ReadSignals(&evo, "signals.par");   
        fflush(stdout); /* useful in MinGW */
    }
    printf("Stop:\n%s\n",  cmaes_TestForTermination(&evo)); /* print termination reason */
    cmaes_WriteToFile(&evo, "all", "allcmaes.dat");         /* write final results */

    /* get best estimator for the optimum, xmean */
    xfinal = cmaes_GetNew(&evo, "xmean"); /* "xbestever" might be used as well */
    cmaes_exit(&evo); /* release memory */ 

    /* do something with final solution and finally release memory */
    free(xfinal); 

    return 0;
}