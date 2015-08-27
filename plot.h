#ifndef PLOT_H_
#define PLOT_H_

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>

class histogram2D {

public:
	histogram2D(double b_x, double b_y);
	int add_to_histogram(double x, double y);	
	int print();
	int print(std::string name);

	double binwidth_x;
	double binwidth_y;

	int dim_X;
	int dim_Y;

	std::vector<std::vector<double>> ARRAY;

private:
	double high_x;
	double low_x;
	double high_y;
	double low_y;

	std::vector<std::vector<double>> print_array;
};

class MMHist {

public:
	MMHist(double b_bins, double b_min, double b_max, int MC_SCALE);
	int wipe_clean();	
	int add_to_histogram(double x);
	int print(double x);
	double efficiency();

private: 
	double binwidth;
	double min;
	double max;
	double bins;
	double events;
	double total_events;
	std::vector<double> histogram;
};


#endif
