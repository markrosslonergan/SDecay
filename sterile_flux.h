#ifndef STERILE_FLUX_H_
#define STERILE_FLUX_H_

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <cstring>


#include "fourmomentum.h"

class initial_sterile {

public: 
	double mass;
	double energy;
	double costhS;
	double phiS;
	fourmomentum labframeP;

	initial_sterile(double M, double E, double in_costhS, double in_phiS);
};

double getEvents(double mS, double mZprime, double events[][2]);

#endif
