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
		print_array.at(floor(ARRAY.at(n).at(0)/binwidth_x) - floor(low_x/binwidth_x)).at(floor(ARRAY.at(n).at(1)/binwidth_y) - floor(low_y/binwidth_y))+=1.0/(ARRAY.size()*binwidth_x*binwidth_y);

	}

	for(n=0;n<dim_X;n++)
	{		
		for(m=0;m<dim_Y;m++)
		{
			// use below for gnuplot splot contour file.
			output_file<<(n+0.5)*binwidth_x<<" "<<(m+0.5)*binwidth_y<<" "<<print_array.at(n).at(m)<<std::endl;	
	
			// use below for gnuplot plot "image" file
			//output_file<<print_array.at(n).at(m)<<" ";	
		}
		output_file<<std::endl;
	}

output_file.close();

return 0;
}


