#include <iostream>
#include <cmath>
#include <vector>

#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "fourmomentum.h" // Defines a class for fourmomenta; makes getting directions and 3-momenta (and
			  // all that) easier.

#include "sterile_flux.h" // Defines getEvents, and describes a class whose objects denote a single 
			  // incoming sterile (mass, fourmomentum).

//#include "decay.h" 	// include functions for turning sterile event objects into
			// observables. Also, defines the structs used to store the variables we care
			// about. I intend this to become a class; and each channel can then become a 
			// derived class.

#include "channel.h"

#define NUMEVENTS 200000

/* ########## Main function ############### */
int main(int argc, char* argv[])
{

double mS = strtof(argv[1],NULL);
double mZprime = strtof(argv[2],NULL); 

static double events[NUMEVENTS][2]; //define the storage for all the events, [0] = E_s, [1] = cos\theta_s
getEvents(mS,mZprime,events); 

const gsl_rng_type * T; // Standard invocation prayer to the RNG
gsl_rng * r;
gsl_rng_env_setup();

T = gsl_rng_default;
r = gsl_rng_alloc (T);

static OBSERVABLES Obs; //This struct is contained in "decay.h"; it specifically gives variables for a two body event (e+,e-)

std::vector<double> model_params; //These are the theoretical parameters which the model needs to know.
model_params.push_back(mZprime);

double phiS = 0.0;


//threebody CHAN(r, model_params);
Zprimeresonance CHAN(r, model_params);

//We enter the main loop over events. For each one, computing the relevant
//observables.
int m; for(m=0;m<=NUMEVENTS-1;m++) 
{

	//The data files I have doesn't provide phi angles to the steriles, so
	//we generate them here. However, I think we want to add the phi to the
	//sterile data. 
	phiS = 2.0*M_PI*gsl_rng_uniform(r);

	//At the moment I package up the sterile parameters into an object here.
	//In future, it would be nice for getEvents to create an array/vector
	//of such objects. Then the loop could just be steping through this
	//array.
	initial_sterile nus(mS, events[m][0], events[m][1], phiS);

	//Here's the paradigm: decay_function(the allocated RNG, pointer to the
	//observables desired, initial sterile, other model parameters e.g.
	//mZprime); 
	//
	// decayfunction(r, &Obs, nus, model_params); 
	//
	// At the moment, you comment out the one you don't want below:
	//
//	resonantZprime_decayfunction(r, &Obs, nus, model_params); 
//	threebody_decayfunction(r, &Obs, nus, model_params); 
	CHAN.decayfunction(nus);
	CHAN.observables(&Obs);

	// The following sterile observables can't be assigned at the channel
	// level anymore... in some sense they are inputs not properties of the
	// outgoing event. We could change this if we want.
	Obs.E_sterile = nus.energy; 	
	Obs.Th_sterile = nus.costhS;

	printf("%.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf\n", Obs.E_sum, Obs.Th_sum, Obs.AngSep, Obs.E_sterile, Obs.Th_sterile, Obs.E_high, Obs.Th_high, Obs.E_low, Obs.Th_low, Obs.FS_AngSep);
}

gsl_rng_free(r);

return 0;
}

