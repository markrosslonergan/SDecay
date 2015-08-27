#include "plot.h"

histogram2D::histogram2D(double b_x, double b_y)
{
	binwidth_x = b_x;
	binwidth_y = b_y;

}

int histogram2D::add_to_histogram(double x, double y)
{

	if(ARRAY.size()==0)
	{ 
		high_x = x; 
		low_x = x; 
		low_y = y; 
		high_y = y;
	}
	else 
	{
		if(high_x<x){high_x=x;}
		else if(low_x>x){low_x = x;}	

		if(high_y<y){high_y=y;}
		else if(low_y>y){low_y=y;}
	}

	std::vector<double> temp; 
	temp.push_back(x);
	temp.push_back(y);

	ARRAY.push_back(temp);

return 0;
}

int histogram2D::print()
{
	dim_X = 1 + floor(high_x/binwidth_x) - floor(low_x/binwidth_x);
	dim_Y = 1 + floor(high_y/binwidth_y) - floor(low_y/binwidth_y);


	std::vector<double> temp;

	int n;
	int m;
	for(m=0;m<dim_Y;m++)
	{		
		temp.push_back(0.0);	
	}
	for(m=0;m<dim_X;m++)
	{		
		print_array.push_back(temp);	
	}

	for(n=0;n<ARRAY.size();n++)
	{
		print_array.at(floor(ARRAY.at(n).at(0)/binwidth_x) - floor(low_x/binwidth_x)).at(floor(ARRAY.at(n).at(1)/binwidth_y) - floor(low_y/binwidth_y))+=1.0/(ARRAY.size()*binwidth_x*binwidth_y);

	}

	for(n=0;n<dim_X;n++)
	{		
		for(m=0;m<dim_Y;m++)
		{
			std::cout<<(n+0.5)*binwidth_x<<" "<<(m+0.5)*binwidth_y<<" "<<print_array.at(n).at(m)<<std::endl;	
		}
		std::cout<<std::endl;
	}



return 0;
}

int histogram2D::print(std::string name)
{

std::ofstream output_file;
output_file.open(name);

	dim_X = 1 + floor(high_x/binwidth_x) - floor(low_x/binwidth_x);
	dim_Y = 1 + floor(high_y/binwidth_y) - floor(low_y/binwidth_y);


	std::vector<double> temp;

	int n;
	int m;
	for(m=0;m<dim_Y;m++)
	{		
		temp.push_back(0.0);	
	}
	for(m=0;m<dim_X;m++)
	{		
		print_array.push_back(temp);	
	}

	for(n=0;n<ARRAY.size();n++)
	{
		//print_array.at(floor(ARRAY.at(n).at(0)/binwidth_x) - floor(low_x/binwidth_x)).at(floor(ARRAY.at(n).at(1)/binwidth_y) - floor(low_y/binwidth_y))+=1.0/(ARRAY.size()*binwidth_x*binwidth_y);
		print_array.at(floor(ARRAY.at(n).at(0)/binwidth_x) - floor(low_x/binwidth_x)).at(floor(ARRAY.at(n).at(1)/binwidth_y) - floor(low_y/binwidth_y))+=1.0;

	}

	for(n=0;n<dim_X;n++)
	{		
		for(m=0;m<dim_Y;m++)
		{
			// use below for gnuplot splot contour file.
			output_file<<(floor(low_x/binwidth_x)+n+0.5)*binwidth_x<<" "<<(floor(low_y/binwidth_y)+m+0.5)*binwidth_y<<" "<<print_array.at(n).at(m)<<std::endl;	
	
			// use below for gnuplot plot "image" file
			//output_file<<print_array.at(n).at(m)<<" ";	
		}
		output_file<<std::endl;
	}

output_file.close();

return 0;
}



MMHist::MMHist(double b_bins, double b_min, double b_max)
{
	events=0.0;
	min= b_min;
	max= b_max;
	bins = b_bins;
	binwidth = (b_max-b_min)/b_bins;	

	int n;
	for(n=0; n<bins; n++)
	{
		histogram.push_back(0.0);
	}

}

int MMHist::add_to_histogram(double x)
{
	int bin_temp = floor(x/binwidth);
	if(bin_temp<histogram.size())
	{
		histogram.at(bin_temp) += 1.0;
	}

	events += 1.0;

return 0;
}

int MMHist::wipe_clean()
{
	int N=0;
	for(N=0;N<histogram.size();N++)
	{
		histogram.at(N) = 0.0;
	}

return 0;
}

int MMHist::print(double x)
{
	int n;
	if(events==0.0)
	{
		for(n=0;n<histogram.size();n++)
		{
			//std::cout<<" "<<0.0;
			std::cout<<x<<" "<<n*binwidth +min<<" "<<0.0<<std::endl;
		}
	}
	else
	{
		for(n=0;n<histogram.size();n++)
		{
			//std::cout<<" "<<histogram.at(n)/events;
			std::cout<<x<<" "<<n*binwidth+min<<" "<<histogram.at(n)/events<<std::endl;
		}
	}
	std::cout<<std::endl;

return 0;	
}
