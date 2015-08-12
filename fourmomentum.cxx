#include "fourmomentum.h"

// Things to do to make this better:
//
// 1. Do some proper error handling in the constructor. Exceptions?  
//
// 2. In the print function include more useful information?  
//
// 3. Overload the constructor to allow other ways to create these objects? For
//    example, we could allow the creation of specific particles with hardcoded
//    masses, by specifying just a 3-momentum, or an energy and a direction... 
//

fourmomentum::fourmomentum(double energy, std::vector<double> momentum)
{
	populate(energy, momentum);
}

//This function allows the declaration of empty fourmometa to be filled by .populate(...) later. USE WITH CARE! Most of the functions would error with such a vector. 
fourmomentum::fourmomentum()
{
	E = 0.0;
	modp=0;
	p.push_back(0.0);
	p.push_back(0.0);
	p.push_back(0.0);
	mass = 0.0;
	
}

int fourmomentum::populate(double energy, std::vector<double> momentum)
{

	E = energy;
	p = momentum;

	//Just check that the input is OK.
	if(p.size() != 3) { std::cout<<"ERROR: 3-momentum wrong size."<<std::endl; }

	//Define the total 3-momentum and the full Minkowski norm (and check it's an on-shell 4-momentum).
	modp = sqrt(p.at(0)*p.at(0) + p.at(1)*p.at(1) + p.at(2)*p.at(2));

	mass = E*E - modp*modp; //Using mass as a temporary variable here. True value created a few lines down.
	if(fabs(mass) < 1e-12){ mass = 0.0; }

	if(mass < 0.0 ){ std::cout<<"ERROR: 4-vector is spacelike. This isn't what we had agreed on!"<<std::endl; }
	else{ mass = sqrt(mass); }


}

int fourmomentum::print(std::string name)
{

std::cout<<"Fourvector '"<<name<<"'"<<" = ("<<E<<", "<<p.at(0)<<", "<<p.at(1)<<", "<<p.at(2)<<"),\t"<<"[Inv. Mass^2: "<<mass<<", Norm of 3-momentum: "<<modp<<"]"<<std::endl; 

return 0;
}

std::vector<double> fourmomentum::direction()
{
	std::vector<double> temp;

	if(modp==0)
	{
		std::cout<<"Cannot compute direction: 3-momentum vanishes."<<std::endl;
		temp.push_back(0.0);
		temp.push_back(0.0);
	}
	else{
		temp.push_back(acos(p.at(2)/modp));
		temp.push_back(atan(p.at(1)/p.at(0)));
	}

return temp;	
}

double fourmomentum::gamma()
{
double temp;
	if(mass==0)
	{ 
		std::cout<<"ERROR: Trying to compute gamma factor for NULL four momentum."<<std::endl; 
		temp = 1e-5;
	}	
	else
	{
		temp = E/mass;
	}

return temp;
}
