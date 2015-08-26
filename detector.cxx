#include "detector.h"


int detector::accept(OBSERVABLES * Obs)
{
	std::cout<<"I'm the accept function of a ghost detector!"<<std::endl;

return -1;
}

muBooNE::muBooNE()
{
	Energy_threshold = 0.1; 	// In GeV
	AngSep_threshold = 30.0; 	// In Degrees
	Energy_ratio_threshold = 0.1; 	// Percentage.
}

int muBooNE::accept(OBSERVABLES * Obs)
{
    if(	Obs->E_sum < Energy_threshold || 
	Obs->FS_AngSep < AngSep_threshold ||
	Obs->E_low/Obs->E_high < Energy_ratio_threshold ) 
    { 
	return REJECTED; 
    }
    else
    {
	return ACCEPTED;
    }

}


int nocuts::accept(OBSERVABLES * Obs)
{
	return 0;
}

