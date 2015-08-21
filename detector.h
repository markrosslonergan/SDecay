#ifndef DETECTOR_H_
#define DETECTOR_H_

#include <iostream>
#include <cmath>
#include <vector>
#include <string>

#include "channel.h"

#define ACCEPTED 0
#define REJECTED 1


class detector {

public: 

	virtual int accept(OBSERVABLES * Obs);

};


class nocuts : public detector {

public:
	int accept(OBSERVABLES * Obs);

};


class muBooNE : public detector {

public:
	muBooNE();
	int accept(OBSERVABLES * Obs);

	double Energy_threshold; //Below this events aren't registered.

	double AngSep_threshold; //Below this electorn positron pairs are seen as a single track.

	double Energy_ratio_threshold; 	// If E_low/E_high is below this the 
					// event is treated as a single track.

};

#endif
