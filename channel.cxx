#include "channel.h"

twoIP_channel::twoIP_channel(gsl_rng * g)
{
	std::vector<double> p;
	p.push_back(0.0);	
	p.push_back(0.0);	
	p.push_back(0.0);	

	IP1.populate(1,p);
	IP2.populate(1,p);
	
	r = g;
}

threebody::threebody(gsl_rng *g, double mZ) : twoIP_channel(g)
{
	model_params.push_back(mZ);
}

Zprimeresonance::Zprimeresonance(gsl_rng *g, double mZ) : twoIP_channel(g)
{
	model_params.push_back(mZ);
}

int threebody::decayfunction(initial_sterile nuS)
{

return 0;
}

int Zprimeresonance::decayfunction(initial_sterile nuS)
{

return 0;
}




