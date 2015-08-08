// THIS CODE PROVIDES THE DECAY FUNCTION FOR 
//
//     nu_s \to Zprime (on-shell) \to eplus eminus
//

#include "decay.h"

// THE FOLLOWING IS FROM 2body/decayer.c.

double rot_boost(double costheta, double phi, double gamma, double FOURVEC[4])
{
	double sintheta = sqrt(1-costheta*costheta);

	double beta = sqrt(1.0-1.0/(gamma*gamma));
	double temp[4];
	temp[0]=FOURVEC[0];
	temp[1]=FOURVEC[1];
	temp[2]=FOURVEC[2];
	temp[3]=FOURVEC[3];

	FOURVEC[0] = gamma*temp[0] + gamma*beta*temp[3];
	FOURVEC[1] = gamma*beta*cos(phi)*sintheta*temp[0] + cos(phi)*costheta*temp[1] - sin(phi)*temp[2] + gamma*cos(phi)*sintheta*temp[3];
	FOURVEC[2] = gamma*beta*sin(phi)*sintheta*temp[0] + sin(phi)*costheta*temp[1] + cos(phi)*temp[2] + gamma*sin(phi)*sintheta*temp[3];
	FOURVEC[3] = gamma*beta*costheta*temp[0] - sintheta*temp[1] + gamma*costheta*temp[3];

return 0.0;
}

double fourvec_costheta(double FOURVEC[4])
{
return FOURVEC[3]/sqrt(FOURVEC[1]*FOURVEC[1]+FOURVEC[2]*FOURVEC[2]+FOURVEC[3]*FOURVEC[3]);
}

double fourvec_cosphi(double FOURVEC[4])
{
	double P = sqrt(FOURVEC[1]*FOURVEC[1]+FOURVEC[2]*FOURVEC[2]+FOURVEC[3]*FOURVEC[3]);
	double cosPhi = FOURVEC[1]/(sqrt(1-fourvec_costheta(FOURVEC)*fourvec_costheta(FOURVEC))*P);
return cosPhi;
}

int resonantZprime_decayfunction(gsl_rng * r, OBSERVABLES * output, initial_sterile nuS, std::vector<double> model_params)
{

double mZprime = model_params.at(0);
double Es = nuS.energy;
double costhS = nuS.costhS;
double phiS = nuS.phiS;
double mS = nuS.mass;

//

double Z_E_srf = (mS*mS+mZprime*mZprime)/(2.0*mS);
double Z_P_srf = sqrt(Z_E_srf*Z_E_srf-mZprime*mZprime);

double Z_phi_srf = 0.0;
double Z_costheta_srf = 1.0;

double S_phi_lf = 0.0;
double S_costheta_lf = 1.0;
double S_E_lf = Es;
double S_gamma = Es/mS;

double Z_FOURVEC[] = {0.0,0.0,0.0,0.0};
double EPLUS_FOURVEC[] = {0.0,0.0,0.0,0.0};
double EMINUS_FOURVEC[] = {0.0,0.0,0.0,0.0};
double SUM_FOURVEC[] = {0.0,0.0,0.0,0.0};

double Z_gamma = 1.0;
double cosalpha, sinalpha;
double cosbeta, sinbeta, beta;
double temp = 0.0;

	// Angles of the Z in the sterile rest frame (srf) are evenly distributed on the sphere.
	Z_phi_srf = 2.0*M_PI*gsl_rng_uniform(r);
	Z_costheta_srf = 2.0*gsl_rng_uniform(r) -1.0;
	
	// The labframe phi and costheta for sterile (S_phi_lf, S_costheta_lf) are taken from input.
	S_phi_lf = phiS;
	S_costheta_lf = costhS;	
	S_E_lf = Es;
	S_gamma = S_E_lf/mS;

	// We define the Z fourvector in the sterile restframe.
	Z_FOURVEC[0] = Z_E_srf;
	Z_FOURVEC[1] = Z_P_srf*sqrt(1.0-Z_costheta_srf*Z_costheta_srf)*cos(Z_phi_srf);
	Z_FOURVEC[2] = Z_P_srf*sqrt(1.0-Z_costheta_srf*Z_costheta_srf)*sin(Z_phi_srf);
	Z_FOURVEC[3] = Z_P_srf*Z_costheta_srf;

	//We boost and rotate to move the Z fourvector into the lab frame.
	rot_boost(S_costheta_lf,S_phi_lf,S_gamma,Z_FOURVEC);
	
//	printf("%.5lf %.5lf %.5lf %.5lf %.5lf\n", Z_FOURVEC[0],  fourvec_costheta(Z_FOURVEC), fourvec_cosphi(Z_FOURVEC), Z_costheta_srf, cos(Z_phi_srf)); 

	Z_gamma = Z_FOURVEC[0]/mZprime;
	
	cosalpha = 2.0*gsl_rng_uniform(r) - 1.0;
	beta = 2.0*M_PI*gsl_rng_uniform(r);
	cosbeta = cos(beta);
	sinalpha = sqrt(1.0-cosalpha*cosalpha);
	sinbeta = sqrt(1.0-cosbeta*cosbeta);

	EPLUS_FOURVEC[0] = mZprime/2.0;
	EPLUS_FOURVEC[1] = (mZprime/2.0)*sinalpha*cosbeta;
	EPLUS_FOURVEC[2] = (mZprime/2.0)*sinalpha*sinbeta;
	EPLUS_FOURVEC[3] = (mZprime/2.0)*cosalpha;

	EMINUS_FOURVEC[0] = mZprime/2.0;
	EMINUS_FOURVEC[1] = -(mZprime/2.0)*sinalpha*cosbeta;
	EMINUS_FOURVEC[2] = -(mZprime/2.0)*sinalpha*sinbeta;
	EMINUS_FOURVEC[3] = -(mZprime/2.0)*cosalpha;

	rot_boost(fourvec_costheta(Z_FOURVEC), acos(fourvec_cosphi(Z_FOURVEC)),Z_gamma, EPLUS_FOURVEC); 
	rot_boost(fourvec_costheta(Z_FOURVEC), acos(fourvec_cosphi(Z_FOURVEC)),Z_gamma, EMINUS_FOURVEC); 

	SUM_FOURVEC[0] = EPLUS_FOURVEC[0] + EMINUS_FOURVEC[0], 
	SUM_FOURVEC[1] = EPLUS_FOURVEC[1] + EMINUS_FOURVEC[1];
	SUM_FOURVEC[2] = EPLUS_FOURVEC[2] + EMINUS_FOURVEC[2];
	SUM_FOURVEC[3] = EPLUS_FOURVEC[3] + EMINUS_FOURVEC[3];


	//OBSERVABLES { double E_sum; double Th_sum; double AngSep; double E_sterile; double Th_sterile; double E_high; double Th_high; double E_low; double Th_low; } OBSERVABLES;
	output->E_sum = SUM_FOURVEC[0];	
	output->Th_sum = acos(fourvec_costheta(SUM_FOURVEC));	
	output->AngSep = acos((EPLUS_FOURVEC[1]*EMINUS_FOURVEC[1] + EPLUS_FOURVEC[2]*EMINUS_FOURVEC[2] + EPLUS_FOURVEC[3]*EMINUS_FOURVEC[3])/( sqrt( EPLUS_FOURVEC[1]*EPLUS_FOURVEC[1] + EPLUS_FOURVEC[2]*EPLUS_FOURVEC[2] + EPLUS_FOURVEC[3]*EPLUS_FOURVEC[3])*sqrt( EMINUS_FOURVEC[1]*EMINUS_FOURVEC[1] + EMINUS_FOURVEC[2]*EMINUS_FOURVEC[2] + EMINUS_FOURVEC[3]*EMINUS_FOURVEC[3])));

	if(EMINUS_FOURVEC[3] > 0 && EPLUS_FOURVEC[3] > 0)
	{
		output->FS_AngSep = (180.0/M_PI)*fabs(atan(EPLUS_FOURVEC[1]/EPLUS_FOURVEC[3]) - atan(EPLUS_FOURVEC[1]/EMINUS_FOURVEC[3]));
	}
	else 
	{
		output->FS_AngSep = 180 - (180.0/M_PI)*fabs(atan(EPLUS_FOURVEC[1]/EPLUS_FOURVEC[3]) - atan(EPLUS_FOURVEC[1]/EMINUS_FOURVEC[3]));
	}

	output->E_sterile = S_E_lf;	
	output->Th_sterile = S_costheta_lf;	
	output->E_high = EPLUS_FOURVEC[0];	
	output->Th_high = fourvec_costheta(EPLUS_FOURVEC);	
	output->E_low = EMINUS_FOURVEC[0];	
	output->Th_low = fourvec_costheta(EMINUS_FOURVEC);	
	if(output->E_high < output->E_low)
	{ 	
		temp = output->E_low; 
		output->E_low = output->E_high; 
		output->E_high = temp; 
		
		temp = output->Th_low;
		output->Th_low = output->Th_high;
		output->Th_high = temp;
	}

//	printf("%.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf\n", output->E_sum, output->Th_sum, output->AngSep, output->E_sterile, output->Th_sterile, output->E_high, output->Th_high, output->E_low, output->Th_low);


return 0;
}

