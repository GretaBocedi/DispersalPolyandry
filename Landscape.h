#pragma once

#include <stdio.h>

class Landscape {
public:
	Landscape();
	~Landscape();
	int pres;
	bool suitable;
	double local_K;
	double *local_K0; //in case of stochasticity
	double eps; //env. stochasticity
	//double *theta; //environmental gradient value  (Atkins & Travis 2010)

private:
};


