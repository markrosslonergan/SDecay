#include "sterile_flux.h"

initial_sterile::initial_sterile(double M, double E, double in_costhS, double in_phiS)
{
	mass = M;
	energy = E;
	costhS = in_costhS; 
	phiS = in_phiS;
	
	double temp[] = {sqrt(energy*energy-M*M)*sqrt(1-in_costhS*in_costhS)*cos(in_phiS), sqrt(energy*energy-M*M)*sqrt(1-in_costhS*in_costhS)*sin(in_phiS), sqrt(energy*energy-M*M)*in_costhS };

	std::vector<double> momentum(temp, temp + sizeof(temp)/sizeof(double));

	labframeP.populate(energy, momentum);
	
}

double getEvents(double mS, double mZprime, double events[][2])
{
	FILE *ptr_file;
    	
	char buf[3000];

	char * pch;

	int n = 1;
	int m = 0;
	char s[100];
	char filename[500] = "../MC/MC_\0";
	sprintf(s,"%.4lf_%.4lf.dat", mS, mZprime);
	strcat(filename,s);
//	printf("Filename: %s\n",filename);
	ptr_file =fopen(filename,"r");

    	if (!ptr_file)
       	{			
		printf("ERROR LOADING MC EVENTS\n");
		exit(1);
	}
    	while (fgets(buf,3000, ptr_file)!=NULL)
	{
		pch = strtok(buf,"\t");
		n=1;
 		while (pch != NULL)
		{
			if(n==3){	//printf("%.7g, ",strtof(pch,NULL));
					events[m][0] = strtof(pch,NULL);	 
				}
			if(n==4){	//printf("%.7g\n",strtof(pch,NULL));
					events[m][1] = strtof(pch,NULL);	 
				}
			pch = strtok(NULL,"\t");
		n++;	
		}
		m++;
	}
	fclose(ptr_file);

//	printf("%s\n",flux_temp);

//printf("Total lines: %d\n", m);

return 0;
}


