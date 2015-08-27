#include <iostream>
#include <cmath>
#include <vector>
#include <unistd.h>

#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "fourmomentum.h" // Defines a class for fourmomenta; makes getting directions and 3-momenta (and
			  // all that) easier.

#include "sterile_flux.h" // Defines getEvents, and describes a class whose objects denote a single 
			  // incoming sterile (mass, fourmomentum).

#include "channel.h"	  // Includes the mother class for two body final state
			  // decays, and the derived classes for the two channels we've studied so far.

#include "plot.h"	  // This includes functions to make 2D histograms.

#include "detector.h"	  // This includes the detector-specific cut functions.

#define NUMEVENTS 200000

int output_distributions(gsl_rng * r, detector * DETECTOR, twoIP_channel * CHAN, double mS, double mZprime);
int migration_matrix(detector * DETECTOR, twoIP_channel * CHAN, double mS);	

/* ########## Main function ############### */
int main(int argc, char* argv[])
{

const gsl_rng_type * T; // Standard invocation prayer to the RNG
gsl_rng * r;
gsl_rng_env_setup();

T = gsl_rng_default;
r = gsl_rng_alloc (T);

double mS = 0.0530; 	 // These are the default values if no command line parameters are defined. 
double mZprime = 0.3600; // This one too.
int channel_flag = 0; 	 // This one too.
int matrix_flag = 0; 	 // This one too.

int c;
opterr = 0;

while ((c = getopt(argc, argv, "m:Z:C:M:")) != -1)
{
switch(c)
{
//
// this is where we create an object for the channel that we want to use.
//
case 'm':
	mS = strtof(optarg,NULL);
	break;
case 'Z':
	mZprime = strtof(optarg,NULL);
	break;
case 'C':
	channel_flag = strtod(optarg,NULL);
	break;
case 'M':
	matrix_flag = strtod(optarg,NULL);
	break;
case '?':
//	std::cout<<"Abandon hope all ye who enter this value. "<<std::endl;
	std::cout<<"Allowed arguments:"<<std::endl;
	std::cout<<"\t-m\tsets the sterile mass. [default = 0.0530]"<<std::endl;
	std::cout<<"\t-Z\tsets the Zprime mass. [default = 0.3600]"<<std::endl;
	std::cout<<"\t-C\tsets the decay channel (0 normal threebody, 1 resonant threebody, 2 generic two body). [default = 0]"<<std::endl;
	std::cout<<"\t-M n\tproduces a migration matrix for oberservable 'n'"<<std::endl;
	return 0;
default:
	std::cout<<"I don't know how you got here."<<std::endl;
	return 0;
}}


//Now we set up the decay channel object.
std::vector<double> model_params; //This should include any theoretical parameters which the model needs to know.

twoIP_channel *CHAN;

if(channel_flag == 1)
{
//	std::cout<<"Resonant threebody."<<std::endl;

	model_params.push_back(mZprime);
	CHAN = new Zprimeresonance(r,model_params); 
}
else if(channel_flag==2)
{
//	std::cout<<"Twobody."<<std::endl;

	model_params.push_back(0.005); //set one decay product massless.
	model_params.push_back(0.010); //set one massive.
	CHAN = new twobody(r,model_params); 
}
else
{
//	std::cout<<"Threebody."<<std::endl;

	model_params.push_back(mZprime);
	CHAN = new threebody(r,model_params); 
}

//We define the detector cuts we would like to use. 
//
detector * DETECTOR;

DETECTOR = new nocuts(); 	// this is a pseudo-detector that just allows every event.
//DETECTOR = new muBooNE(); 	// this is microboone.


if(matrix_flag == 0)
{
	//This is the meat of the old program
	output_distributions(r, DETECTOR, CHAN, mS, mZprime);
}
else
{
	migration_matrix(DETECTOR, CHAN, mS);	
}


//Cleaning up.
gsl_rng_free(r);
delete CHAN;

return 0;
}


int migration_matrix(detector * DETECTOR, twoIP_channel * CHAN, double mS)
{

//For now we assume that all steriles are on axis. cos=1 phi=0
double Emin=0.05;
double Emax=10.0;
double number_of_bins = 40;

static OBSERVABLES Obs; //This struct is contained in "decay.h"; it specifically gives variables for a two body event (e+,e-)

//These should never get changed and never used, but I thought it safe to fill them with something.
Obs.Th_sterile = 0.0;
Obs.E_sterile = 0.0;

int m; 
int MC_SCALE = 40000;

MMHist EsumHist(100.0,0.0,5.0);

double Es = 1.0;
for(Es=Emin; Es<Emax+1e-5; Es+=(Emax-Emin)/number_of_bins)
{

if(Es>mS)
{
	initial_sterile nus(mS, Es, 1.0-1e-7, 0.0); //I'm a little wary of putting the steriles exactly on axis... why? Test it?

	EsumHist.wipe_clean();

	for(m=0;m<MC_SCALE;m++)
	{
		//We call the appropriate functions from the channels.
		CHAN->decayfunction(nus);
		CHAN->observables(&Obs);

		if(DETECTOR->accept(&Obs)==ACCEPTED)
		{
			EsumHist.add_to_histogram(Obs.E_sum);
		}
	}
}
	EsumHist.print(Es); //this prints zeros if no events have been added.
}	

return 0;
}



int output_distributions(gsl_rng * r, detector * DETECTOR, twoIP_channel * CHAN, double mS, double mZprime)
{

static double events[NUMEVENTS][2]; //define the storage for all the events, [0] = E_s, [1] = cos\theta_s
getEvents(mS,mZprime,events); 

//We can also make a histogram suitable for gnuplot.  The two arguments are the
//binwidths in the x and y variables (in this case Esum and the foreshortened
//angular separation).
histogram2D HIST_ESUM_FSANGSEP(0.01,1.0);

//and one for the total desposited energy against the high/low energy ratio.
histogram2D HIST_ESUM_EHIGHLOWRATIO(0.01,0.01);

//For angular separation against foreshortened angular separation.
histogram2D HIST_ANGSEP_FSANGSEP(1.0,1.0);

//We enter the main loop over events. For each one, computing the relevant
//observables.
double phiS = 0.0;
static OBSERVABLES Obs; //This struct is contained in "decay.h"; it specifically gives variables for a two body event (e+,e-)
int m; for(m=0;m<=NUMEVENTS-1;m++) 
{
	//The data files I have don't provide phi angles for the steriles, so
	//we generate them here. However, I think we want to add the phi to the
	//sterile data. 
	phiS = 2.0*M_PI*gsl_rng_uniform(r);

	//At the moment I package up the sterile parameters into an object here.
	//In future, it would be nice for getEvents to create an array/vector
	//of such objects. Then the loop could just be steping through this
	//array.
	initial_sterile nus(mS, events[m][0], events[m][1], phiS);

	//We call the appropriate functions from the channels.
	CHAN->decayfunction(nus);
	CHAN->observables(&Obs);

	// The following sterile observables can't be assigned at the channel
	// level anymore... in some sense they are inputs, not properties of the
	// outgoing event, so I'm not sure I think this is a problem.
	Obs.E_sterile = nus.energy; 	
	Obs.Th_sterile = nus.costhS;

	if(DETECTOR->accept(&Obs)==ACCEPTED)
	{
		printf("%.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf\n", Obs.E_sum, Obs.Th_sum, Obs.AngSep, Obs.E_sterile, Obs.Th_sterile, Obs.E_high, Obs.Th_high, Obs.E_low, Obs.Th_low, Obs.FS_AngSep);
	
		HIST_ESUM_FSANGSEP.add_to_histogram(Obs.E_sum,Obs.FS_AngSep);
		HIST_ESUM_EHIGHLOWRATIO.add_to_histogram(Obs.E_sum,Obs.E_low/Obs.E_high);
		HIST_ANGSEP_FSANGSEP.add_to_histogram(Obs.AngSep,Obs.FS_AngSep);
	}
	
}

	HIST_ESUM_FSANGSEP.print("data/Esum_FSangularsep.dat");
	HIST_ESUM_EHIGHLOWRATIO.print("data/Esum_EnergyRatio.dat");
	HIST_ANGSEP_FSANGSEP.print("data/Angularsep_FSangularsep.dat");

return 0;
}

