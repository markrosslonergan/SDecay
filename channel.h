#ifndef CHANNEL_H_
#define CHANNEL_H_

#include <iostream>
#include <cmath>
#include <vector>

#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "fourmomentum.h" // defines class fourmomentum
#include "sterile_flux.h" // defines class initial_sterile

class twoIP_channel { //This is the mother class for all decay channels (into two Ionising Particles)

public:
	twoIP_channel(gsl_rng * g);

	fourmomentum IP1;	//first outgoing particle 4 momentum.
	fourmomentum IP2;	//second outgoing particle 4 momentum.

	gsl_rng * r;
	std::vector<double> model_params;
};


/* ###############################
   
   Below here we have a derived class for each channel

   ############################### */


//This is the nu_s \to \nu e+ e- channel (off-shell Zprime).
class threebody : public twoIP_channel {

public:
	threebody(gsl_rng * g, double mass);
	int decayfunction(initial_sterile nuS);

}; 

//This is the nu_s \to \nu Zprime \to \nu e+ e- channel (on-shell Zprime).
class Zprimeresonance : public twoIP_channel {

public: 
	Zprimeresonance(gsl_rng * g, double mass);
	int decayfunction(initial_sterile nuS);

}; 

#endif
