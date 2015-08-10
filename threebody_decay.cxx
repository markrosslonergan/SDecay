//THE FOLLOWING IS FROM THE DECAY DECAYER (Oh, what? That sounds redundant, eh? FUCK YOU.)
//
// The code computes the normal three body decay:
// 
//		nu_s \to \nu eplus eminus 
//
// mediated by the off-shell Zprime.

#include "decay.h"

double pdf_function(double x, double y, double mS, double mZprime, void * pointer)
{
	double mu_s  = mS/mZprime;
	double alpha = mu_s*mu_s/(1.0-mu_s*mu_s);

	double invnorm;

	if(alpha < 0.01)
	{
		double invnorm_perturb_0 = (1.0/(1.0+alpha))*(7.0/4.0 + 41.0/60.0*alpha); 

		double invnorm_perturb_rest = -0.18333*pow(alpha,2.0)+0.22857*pow(alpha,3.0)-0.23274*pow(alpha,4.0)+0.22421*pow(alpha,5.0)-0.21190*pow(alpha,6.0)+0.19899*pow(alpha,7.0)-0.18662*pow(alpha,8.0)+0.17517*pow(alpha,9.0)-0.16475*pow(alpha,10.0)+0.15531*pow(alpha,11.0);
		invnorm = invnorm_perturb_0+invnorm_perturb_rest;
	}
	else 
	{		
		invnorm = (3.0/(2.0*pow(alpha,3.0)))*(2+alpha-3*pow(alpha,2.0))/(1.0+alpha) + (4.0*alpha*alpha-3.0)*(log(1+alpha)/log(exp(1.0)))/pow(alpha,4.0);
	}

	double ret = (1.0/invnorm)*x*(4-x*x)/((1.0+alpha*x)*(1.0+alpha*x));

//	printf("inv. norm: %.5lf\talpha: %.5lf\tret: %.5lf\n",invnorm,alpha,ret);

	if(ret<0){ ret = 0.0; }

return ret; 
}

struct PDF_CHOICE choose_from_pdf(gsl_rng * r,double mS, double mZprime, double (*pdf)(double, double, double, double, void *))
{
	double mu_s  = mS/mZprime;
	double alpha = mu_s*mu_s/(1-mu_s*mu_s);

//	double PDF_MAX = 3.0/2.0; //pdf_function_test
	double PDF_MAX = 1.8; //pdf_function

	double x = gsl_rng_uniform(r);
	double y = -1.0 + 2.0*gsl_rng_uniform(r);
	double phi = 2*M_PI*gsl_rng_uniform(r);
	double z = (PDF_MAX+0.01)*gsl_rng_uniform(r);

	while(pdf(x,y,mS,mZprime,NULL)<=z)
	{
//		printf("I tried!\n");
		x = gsl_rng_uniform(r);
		y = -1.0 + 2.0*gsl_rng_uniform(r);
		z = (PDF_MAX+0.01)*gsl_rng_uniform(r);
	}

//	printf("I succeeded!\n");

	//printf("%.5lf %.5lf %.5lf %.5lf\n",x,y,z,pdf(x,y,mS,mZprime,NULL));
	struct PDF_CHOICE package;
	package.Enu = mS*x/2.0;
	package.cosThnu = y;
	package.Phinu = phi;

return package;
}

int rotor(double theta, double phi, double vector[3])
{
	double x=vector[0];
	double y=vector[1];
	double z=vector[2];

	double rdotn1 = x*cos(phi) + y*sin(phi);
	double rdotn2 = z;
	double rdotn3 = -x*sin(phi) + y*cos(phi);

	vector[0]=(cos(theta)*rdotn1 + sin(theta)*rdotn2)*cos(phi) - rdotn3*sin(phi);
	vector[1]=(cos(theta)*rdotn1 + sin(theta)*rdotn2)*sin(phi) + rdotn3*cos(phi);
	vector[2]=-sin(theta)*rdotn1 + cos(theta)*rdotn2;

return 0;
}

int drawRestFrameDist(gsl_rng * r, double mS, double mZprime, double output[3])
{
	struct PDF_CHOICE choice=choose_from_pdf(r,mS,mZprime,pdf_function);
	output[0]=choice.Enu;
	output[1]=choice.cosThnu;
	output[2]=choice.Phinu;
return 0;
}

int computeLabFrameVariables(OBSERVABLES * output, double mS, double Es, double costhS, double phiS, double restFrameParams[3])
{
	double Enu = restFrameParams[0];
	double Th = acos(restFrameParams[1]);
	double Phi = restFrameParams[2];

	double me = 0.0;// THIS CAUSES ERRORS! 5.11e-6; //GeV
	
	double Ee = (mS - Enu)/2.0;
	double Pe = sqrt(Ee*Ee-me*me);
	double beta = sqrt(1-mS*mS/(Es*Es));
	double gamma = 1.0/sqrt(1.0-beta*beta);

//printf("%.5lf %.5lf %.5lf %.5lf %.5lf %.5lf\n", Enu, Ee, Pe, beta, gamma, mS); 

	double alpha = 2.0*acos( Enu/(2.0*Pe) );
	double theta_plus = M_PI - Th - alpha/2.0;
	double theta_minus = M_PI - Th + alpha/2.0;

	double Pplus_E = gamma*(Ee + beta*Pe*cos(theta_plus));
	double Pminus_E = gamma*(Ee + beta*Pe*cos(theta_minus));
	double Pplus_x = Pe*sin(theta_plus)*cos(Phi);
	double Pminus_x = Pe*sin(theta_minus)*cos(Phi);
	double Pplus_y = Pe*sin(theta_plus)*sin(Phi);
	double Pminus_y = Pe*sin(theta_minus)*sin(Phi);
	double Pplus_z = gamma*(Pe*cos(theta_plus) + beta*Ee);
	double Pminus_z = gamma*(Pe*cos(theta_minus) + beta*Ee);

	double Pee[] = {(Pplus_x + Pminus_x)/2.0, (Pplus_y + Pminus_y)/2.0, (Pplus_z + Pminus_z)/2.0};
	double Pplus[] = {Pplus_x, Pplus_y, Pplus_z};
	double Pminus[] = {Pminus_x, Pminus_y, Pminus_z};
	rotor(acos(costhS),phiS,Pee);
	rotor(acos(costhS),phiS,Pplus);
	rotor(acos(costhS),phiS,Pminus);


	output->E_sum = Pplus_E + Pminus_E; // energy is unaffected by rotation.
	output->Th_sum = (180.0/M_PI)*acos(Pee[2]/sqrt(Pee[0]*Pee[0] + Pee[1]*Pee[1] + Pee[2]*Pee[2] )); 
	output->AngSep = (180.0/M_PI)*acos((Pplus_x*Pminus_x + Pplus_y*Pminus_y + Pplus_z*Pminus_z)/(sqrt(Pplus_x*Pplus_x + Pplus_y*Pplus_y + Pplus_z*Pplus_z)*sqrt(Pminus_x*Pminus_x + Pminus_y*Pminus_y + Pminus_z*Pminus_z))); // opening angle is unaffected by rotation

	if(Pminus_z > 0 && Pplus_z > 0)
	{
		output->FS_AngSep = (180.0/M_PI)*fabs(atan(Pplus_x/Pplus_z) - atan(Pminus_x/Pminus_z));
	}
	else 
	{
		output->FS_AngSep = 180- (180.0/M_PI)*fabs(atan(Pplus_x/Pplus_z) - atan(Pminus_x/Pminus_z));
	}

	if(Pplus_E < Pminus_E)
	{
		output->E_low = Pplus_E;  // E_low
		output->E_high = Pminus_E; // E_high

		output->Th_low = (180.0/M_PI)*acos(Pplus[2]/sqrt(Pplus[0]*Pplus[0] +Pplus[1]*Pplus[1] +Pplus[2]*Pplus[2]));
		output->Th_high = (180.0/M_PI)*acos(Pminus[2]/sqrt(Pminus[0]*Pminus[0] +Pminus[1]*Pminus[1] +Pminus[2]*Pminus[2]));
	}
	else 
	{
		output->E_low = Pminus_E; // E_low
		output->E_high = Pplus_E;  // E_high

		output->Th_low = (180.0/M_PI)*acos(Pminus[2]/sqrt(Pminus[0]*Pminus[0] +Pminus[1]*Pminus[1] +Pminus[2]*Pminus[2]));
		output->Th_high = (180.0/M_PI)*acos(Pplus[2]/sqrt(Pplus[0]*Pplus[0] +Pplus[1]*Pplus[1] +Pplus[2]*Pplus[2]));
	}

return 0;
}

// THIS IS JUST A WRAPPER FUNCTION TO GET ALL THE PREVIOUS CODE IN A MODULAR FORMAT.
//
int threebody_decayfunction(gsl_rng *r, OBSERVABLES * output, initial_sterile nuS, std::vector<double> model_params)
{
	double mZprime = model_params.at(0);
	double restFrameParams[] = {0.0,0.0,0.0};
	drawRestFrameDist(r,nuS.mass,mZprime,restFrameParams); //this will populate the doubles.
	computeLabFrameVariables(output, nuS.mass, nuS.energy, nuS.costhS, nuS.phiS, restFrameParams);
	output->E_sterile = nuS.energy;
	output->Th_sterile = nuS.costhS;

return 0;
}



