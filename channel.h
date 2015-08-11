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

typedef struct OBSERVABLES { 
	double E_sum; 	
	double Th_sum; 
	double AngSep; 
	double E_sterile; 
	double Th_sterile; 
	double E_high; 
	double Th_high; 
	double E_low; 
	double Th_low; 
	double FS_AngSep; //The foreshortened angular separation.
	} OBSERVABLES;

class twoIP_channel { //This is the mother class for all decay channels (into two Ionising Particles)

public:
	twoIP_channel(gsl_rng * g, std::vector<double> input);

	fourmomentum IP1;	//first outgoing particle 4 momentum.
	fourmomentum IP2;	//second outgoing particle 4 momentum.

	gsl_rng * r;
	std::vector<double> model_params;

	int observables(OBSERVABLES * output);
};


/* ###############################
   
   Below here we have a derived class for each channel

   ############################### */


//This is the nu_s \to \nu e+ e- channel (off-shell Zprime).
class threebody : public twoIP_channel {

public:
	threebody(gsl_rng * g, std::vector<double> input);
	int decayfunction(initial_sterile nuS);

	struct PDF_CHOICE { 
		double Enu; 
		double cosThnu; 
		double Phinu; };

//	typedef double (*threebody_pdf_function)(double, double, double, double, void *);

private:
	int computeLabFrameVariables(double mS, double Es, double costhS, double phiS, double restFrameParams[3]);
	double pdf_function(double x, double y, double mS, double mZprime, void * pointer);
//	struct PDF_CHOICE choose_from_pdf(gsl_rng * r, double mS, double mZprime, threebody_pdf_function pdf);
	int rotor(double theta, double phi, double vector[3]);
	int drawRestFrameDist(gsl_rng * r, double mS, double mZprime, double output[3]);
}; 

//This is the nu_s \to \nu Zprime \to \nu e+ e- channel (on-shell Zprime).
class Zprimeresonance : public twoIP_channel {

public: 
	Zprimeresonance(gsl_rng * g, std::vector<double> input);
	int decayfunction(initial_sterile nuS);

private:
	double fourvec_costheta(double FOURVEC[4]);
	double fourvec_cosphi(double FOURVEC[4]);
	double rot_boost(double costh, double phi, double gam, double FOURVEC[4]);	
}; 

#endif
