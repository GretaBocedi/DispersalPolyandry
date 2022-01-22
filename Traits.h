#pragma once

#include <stdio.h>

#include "Parameters.h"

class Traits
{
public:
	Traits();
	~Traits();

#if ONECHROM
#else
	double *a_mat; //re-mating slope (probability of mating = exp(-a_mat * n_mates)
	double *ep_f, *ep_m; //female and male emigration probability 
	double *dist_f, *dist_m; //female and male mean dispersal distance
#endif	

	//genotypic values
	double *g_a_mat;
	double *g_ep_f, *g_ep_m;
	double *g_dist_f, *g_dist_m;

	//phenotypic values
	double *p_a_mat;
	double *p_ep_f, *p_ep_m;
	double *p_dist_f, *p_dist_m;

	void initialise(Parameters);
	void deleteTraits(void);

private:
};


