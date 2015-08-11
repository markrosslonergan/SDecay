#include "channel.h"

twoIP_channel::twoIP_channel(gsl_rng * g, std::vector<double> input_params)
{
	std::vector<double> p;
	p.push_back(0.0);	
	p.push_back(0.0);	
	p.push_back(0.0);	

	IP1.populate(1,p);
	IP2.populate(1,p);
	
	model_params = input_params;

	r = g;
}

int twoIP_channel::observables(OBSERVABLES * output)
{
	//OBSERVABLES { double E_sum; double Th_sum; double AngSep; double E_sterile; double Th_sterile; double E_high; double Th_high; double E_low; double Th_low; double FS_AngSep; } OBSERVABLES;

	fourmomentum sum;
	std::vector<double> sum_p;
	sum_p.push_back(IP1.p.at(0) + IP2.p.at(0));
	sum_p.push_back(IP1.p.at(1) + IP2.p.at(1));
	sum_p.push_back(IP1.p.at(2) + IP2.p.at(2));
	
	sum.populate(IP1.E + IP2.E, sum_p);

	output->E_sum = sum.E;	
	output->Th_sum =(180.0/M_PI)*sum.direction().at(0);	

	output->AngSep =(180.0/M_PI)*acos((IP1.p.at(0)*IP2.p.at(0) + IP1.p.at(1)*IP2.p.at(1) + IP1.p.at(2)*IP2.p.at(2))/(IP1.modp*IP2.modp)); 

//	IP1.print("IP1");
//	IP2.print("IP2");

	if(IP1.p.at(2) > 0 && IP2.p.at(2) > 0)
	{
		output->FS_AngSep = (180.0/M_PI)*fabs(atan(IP1.p.at(0)/IP1.p.at(2)) - atan(IP2.p.at(0)/IP2.p.at(2)));
	}
	else 
	{
		output->FS_AngSep = 180 - (180.0/M_PI)*fabs(atan(IP1.p.at(0)/IP1.p.at(2)) - atan(IP2.p.at(0)/IP2.p.at(2)));
	}

	output->E_high = IP1.E;	
	output->Th_high = (180.0/M_PI)*IP1.direction().at(0);	
	output->E_low = IP2.E;	
	output->Th_low = (180.0/M_PI)*IP2.direction().at(0);	

	double temp = 0.0;
	if(output->E_high < output->E_low)
	{ 	
		temp = output->E_low; 
		output->E_low = output->E_high; 
		output->E_high = temp; 
		
		temp = output->Th_low;
		output->Th_low = output->Th_high;
		output->Th_high = temp;
	}

//	IP1.print("IP1");
//	IP2.print("IP2");

return 0;
}


/* ###############################
   
   Below here we have a section for each channel.

   ############################### */


//This is the nu_s \to \nu e+ e- channel (off-shell Zprime).

threebody::threebody(gsl_rng *g, std::vector<double> input) : twoIP_channel(g, input)
{
	if(model_params.size() != 1)
	{ 
		std::cout<<"ERROR: threebody decay channel set up with too many parameters."<<std::endl; 
	}
}

int threebody::decayfunction(initial_sterile nuS)
{
	double mZprime = model_params.at(0);
	double restFrameParams[] = {0.0,0.0,0.0};
	drawRestFrameDist(r,nuS.mass,mZprime,restFrameParams); //this will populate the doubles.
	computeLabFrameVariables(nuS.mass, nuS.energy, nuS.costhS, nuS.phiS, restFrameParams);
return 0;
}

int threebody::computeLabFrameVariables(double mS, double Es, double costhS, double phiS, double restFrameParams[3])
{
	double Enu = restFrameParams[0];
	double Th = acos(restFrameParams[1]);
	double Phi = restFrameParams[2];

//	std::cout<<"Enu: "<<Enu<<" cosThnu: "<<cos(Th)<<" PhiNu: "<<Phi<<std::endl;
//	std::cout<<"Es: "<<Es<<" cosThS: "<<costhS<<" phiS: "<<phiS<<std::endl;

	double me = 0.0;// THIS CAUSES ERRORS! 5.11e-6; //GeV
	
	double Ee = (mS - Enu)/2.0;
	double Pe = sqrt(Ee*Ee-me*me);
	double beta = sqrt(1-mS*mS/(Es*Es));
	double gamma = 1.0/sqrt(1.0-beta*beta);

//	printf("Enu: %.5lf Ee: %.5lf Pe: %.5lf beta: %.5lf gamma: %.5lf mS: %.5lf\n", Enu, Ee, Pe, beta, gamma, mS); 

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

//	std::cout<<"PplusNorm: "<<Pplus_E*Pplus_E - Pplus_x*Pplus_x - Pplus_y*Pplus_y - Pplus_z*Pplus_z<<std::endl;
//	std::cout<<"PmiusNorm: "<<Pminus_E*Pminus_E - Pminus_x*Pminus_x - Pminus_y*Pminus_y - Pminus_z*Pminus_z<<std::endl;

//	std::cout<<"Pp1: "<<Pplus[0]<<" Pp2: "<<Pplus[1]<<" Pp3: "<<Pplus[2]<<std::endl;
//	std::cout<<"Pm1: "<<Pminus[0]<<" Pm2: "<<Pminus[1]<<" Pm3: "<<Pminus[2]<<std::endl;

//	std::vector<double> Vec_pplus(Pplus, Pplus + sizeof(Pplus)/sizeof(double));
//	std::vector<double> Vec_pminus(Pminus, Pminus + sizeof(Pminus)/sizeof(double));

//	IP1.populate(Pplus_E, Vec_pplus);
//	IP1.print("pre-rot IP1");
//	IP2.populate(Pminus_E, Vec_pminus);
//	IP2.print("pre-rot IP2");

	rotor(acos(costhS),phiS,Pee);
	rotor(acos(costhS),phiS,Pplus);
	rotor(acos(costhS),phiS,Pminus);
//	std::cout<<"Rotated!..."<<std::endl;
//	std::cout<<"Pp1: "<<Pplus[0]<<" Pp2: "<<Pplus[1]<<" Pp3: "<<Pplus[2]<<std::endl;
//	std::cout<<"Pm1: "<<Pminus[0]<<" Pm2: "<<Pminus[1]<<" Pm3: "<<Pminus[2]<<std::endl;

	std::vector<double> Vec_pplus2(Pplus, Pplus + sizeof(Pplus)/sizeof(double));
	std::vector<double> Vec_pminus2(Pminus, Pminus + sizeof(Pminus)/sizeof(double));

	IP1.populate(Pplus_E, Vec_pplus2);
//	IP1.print("post-rot IP1");
	IP2.populate(Pminus_E, Vec_pminus2);
//	IP2.print("post-rot IP2");

//std::cout<<std::endl;

return 0;
}

double threebody::pdf_function(double x, double y, double mS, double mZprime, void * pointer)
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


int threebody::rotor(double theta, double phi, double vector[3])
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

int threebody::drawRestFrameDist(gsl_rng * r, double mS, double mZprime, double output[3])
{

	double mu_s  = mS/mZprime;
	double alpha = mu_s*mu_s/(1-mu_s*mu_s);

//	double PDF_MAX = 3.0/2.0; //pdf_function_test
	double PDF_MAX = 1.8; //pdf_function

	double x = gsl_rng_uniform(r);
	double y = -1.0 + 2.0*gsl_rng_uniform(r);
	double phi = 2*M_PI*gsl_rng_uniform(r);
	double z = (PDF_MAX+0.01)*gsl_rng_uniform(r);

	while(threebody::pdf_function(x,y,mS,mZprime,NULL)<=z)
	{
//		printf("I tried!\n");
		x = gsl_rng_uniform(r);
		y = -1.0 + 2.0*gsl_rng_uniform(r);
		z = (PDF_MAX+0.01)*gsl_rng_uniform(r);
	}

//	printf("I succeeded!\n");

	//printf("%.5lf %.5lf %.5lf %.5lf\n",x,y,z,pdf(x,y,mS,mZprime,NULL));
	struct threebody::PDF_CHOICE choice;
	choice.Enu = mS*x/2.0;
	choice.cosThnu = y;
	choice.Phinu = phi;

	output[0]=choice.Enu;
	output[1]=choice.cosThnu;
	output[2]=choice.Phinu;

return 0;
}

//This is the nu_s \to \nu Zprime \to \nu e+ e- channel (on-shell Zprime).

Zprimeresonance::Zprimeresonance(gsl_rng *g, std::vector<double> input) : twoIP_channel(g, input)
{
	if(model_params.size() != 1)
	{
	 	std::cout<<"ERROR: Zprime resonance channel set up with too many parameters."<<std::endl;
	}
}

int Zprimeresonance::decayfunction(initial_sterile nuS)
{

return 0;
}




