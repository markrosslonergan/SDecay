#ifndef DECAY_H_
#define DECAY_H_

#include <iostream>
#include <cmath>
#include <vector>

#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "sterile_flux.h"

typedef struct OBSERVABLES { 
	double E_sum; 	
	double Th_sum; 
	double AngSep; 
	double E_sterile; 
	double Th_sterile; 
	double E_high; 
	double Th_high; 
	double E_low; 
	double Th_low; } OBSERVABLES;

struct PDF_CHOICE { 
	double Enu; 
	double cosThnu; 
	double Phinu; };

// THE FOLLOWING PART IS FROM THE ON-SHELL ZPRIME CODE

double rot_boost(double costheta, double phi, double gamma, double FOURVEC[4]);
double fourvec_costheta(double FOURVEC[4]);
double fourvec_cosphi(double FOURVEC[4]);

int resonantZprime_decayfunction(gsl_rng * r, OBSERVABLES * output, initial_sterile nuS, std::vector<double> model_params);


// THE STUFF BELOW IS FROM THE ORIGINAL MINIBOONE CODE (threebody: eplus eminus neutrino).

double pdf_function(double x, double y, double mS, double mZprime, void * pointer);
struct PDF_CHOICE choose_from_pdf(gsl_rng * r, double mS, double mZprime, double (*pdf)(double, double, double, double, void *));
int rotor(double theta, double phi, double vector[3]);
int drawRestFrameVariable(gsl_rng * r, double mS, double mZprime, double output[3]);
int computeLabFrameVariables(OBSERVABLES * output, double mS, double Es, double costhS, double phiS, double restFrameParams[3]);

int threebody_decayfunction(gsl_rng * r, OBSERVABLES * output, initial_sterile nuS, std::vector<double> model_params);


#endif
