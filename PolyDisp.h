#pragma once

#define CLUSTER 0


#include <stdio.h>
#include <stdlib.h>
#if CLUSTER 
#include <unistd.h>
#else
#include <tchar.h> 
#include <direct.h>
#include <io.h>
#endif
#include <vector>
#include <set>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <math.h>
#include <numeric>
#include <time.h>
#include <random>
#include <iterator>

#include "Parameters.h"
#include "Landscape.h"
#include "Population.h"


using namespace std;

//-----------------------------------------------------------------------------
const double PI = 3.141592654;

clock_t extime;

Parameters para;
string dir, dirOut;
int r, g; // counters for replicates and generations
int cy_max; //max y_colonised 
int cy_min; //minimum y
int Ntot;
bool mat_trials;
ofstream par, inds, offs, pops, trait, muts,  ranInd, popmut;

#if STOPEVO
int Npop;
double mean_disp;
#endif


//-----------------------------------------------------------------------------
Landscape ***land;
Population ***pop;

vector<Individuals> random_offs;

//-----------------------------------------------------------------------------
// Random numbers generators
// seed random number generator
std::random_device rd;
std::mt19937 rdgen(rd());

//-----------------------------------------------------------------------------
std::uniform_int_distribution<> uniint(0, 1);
std::uniform_real_distribution<> unireal(0.0, 1.0);
std::bernoulli_distribution Bern(0.5); 

//probability of population extinciton
std::bernoulli_distribution extinct(para.p_cell_ext);
//temporal environmental stochasticity
std::normal_distribution<> t_stoch(0.0, para.estd);

//Sample locus for adaptive trait mutations
std::uniform_int_distribution<> uni_loci(0, para.L * 2 - 1);
//Recombiantion probability for adaptive traits
std::bernoulli_distribution recomb(para.rec_probability);

//Recombination position increment for neutral alleles
std::geometric_distribution<> geom(para.rec_probability);

//Number of neutral mutations
std::poisson_distribution<> n_neutmut(2.0 * (double)para.nL * para.mu_neutral);
//neutral mutation position
std::uniform_int_distribution<> neutr_position(0, 2 * para.nL - 1);
//Distribution for initialising neutral markers and for their mutations
std::uniform_real_distribution<> neutral(-1000.0, 1000.0);

//Distributions for sampling deleterious mutations
//sample no. of crossovers
std::poisson_distribution<> crossn(para.R);
//Deleterious mutation position (and crossovers positions)
std::uniform_real_distribution<> position(0.0, para.R);
//Selection coefficient for mildly deleterious mutations
std::gamma_distribution<> s_mild(1.0, para.mean_sd); //(shape, mean / shape)
//Dominance coefficient for mildly deleterious mutations 
double k = -log(2.0 * para.mean_hd) / para.mean_sd;

//adaptive mutations (per individual/ per trait)
std::poisson_distribution<> n_amut(2.0 * (double)para.L * para.mu);

//Number of offspring
std::poisson_distribution<> pois(para.fec);


//-----------------------------------------------------------------------------
// Functions declaration
const string Int2Str(const int x);
const string Float2Str(const double x);
void RunModel(void);
void landscape(int, int, double);
void initialisation(int, int, int, int, int, int);
void mating(void);
void birth(void);
void env_stochasticity(int xx, int yy);
void local_extinction(int xx, int yy);
//void generateIDdata(void);
void betweenPop_mating(void);
void inheritance(Individuals*, Individuals, Individuals);
void adaptiveMutations(Individuals*, std::normal_distribution<>, std::normal_distribution<>, std::normal_distribution<>, std::normal_distribution<>);
void premating_dispersal(void);
void postmating_dispersal(void);
void survival(void);
void outPop_header(void);
void outTrait_header(void);
void outMut_header(void);
void outInds_header(void);
void outOffs_header(void);
void outRanInds_header(void);
void outPopMut_header(void);