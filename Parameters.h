#pragma once

#include <fstream>

#define ONECHROM 0 //if 1, deteleterious mutations and trait alleles are on the same chromosome
#define MUTATIONLOAD 1
#define STOPEVO 0
#define POLY 1 //whether female multiple mating is considered or not

using namespace std;

class Parameters
{
public:
	Parameters();
	~Parameters();
	//Simulation
	int SimNr;
	int rep; //replicates
	int gen; //generations
	int out_int; //generations interval for population and trait outputs 
	int out_start;
	int ind_interval; //generations interval for individual outputs 
	int indRand_interval; //generation interval for random mating and output for ID
	int PopMut_interval;
	bool out_mutations;

	//Landscape
	bool spat_heterog; //spatial heterogeneity in K U[Kmin, Kmax]
	int resol; //resolution (m)
	int x_max;
	int y_max;
	int cutoff; //distance from the front after which rows should be cut off
	double prop_suitable; //proportion of suitable cells 
	double K, Kmin, Kmax; //carrying capacity (individuals/ha)

	int exp_start;
	int start_dispEvol;
	int max_age; //max. age of a population to be considered part of the front

	//temporal environmental stochasticity
	bool env_stoch;
	double estd, ac; //amplitude and autocorrelation for time series (Kmin will be zero; K needs to be set to mean K)
	double p_cell_ext; //per-cell extiction probability
	
	////Environmental gradient
	//bool grad_shifts;
	//bool niche; //indicates if species has a niche or not
	//int shift_start, shift_end;
	//double grad_inc;
	//double shift_rate;
	//double opt; // environmental optimum
	//double width; //niche width

	//Initialisation
	int min_seedX;
	int min_seedY;
	int max_seedX;
	int max_seedY;

	//Traits
	bool sexDisp; //TRUE = sex-limited dispersal
	bool nearest; //TRUE = neartest neighbour dispersal 
	int L; //number of loci for each trait
	int dispEvol; //0 = no dispersal evolution; 1 = emigration p. evolving; 2 = disp. distance evolving; 3 = both evolving
	int polyEvol; //0 = no polyandry evolution; 1 = polyandry evolving
	double *genot_mean; // initial genotypic mean for each trait (a, ep_f, ep_m, dist_f, dist_m)
	double *genot_std; // initial genotypic standard deviation for each trait (a, ep_f, ep_m, dist_f, dist_m)
	double *genot_std2; //gentotypic standard deviation for when traits start evolving - DETERMINE MUTATION SD
	double mu; // haploid per allele mutation rate for each trait
	double mu_neutral;
	double *mu_std; // standard deviation for mutational effects for each trait (a, ep_f, ep_m, dist_f, dist_m)
	double rec_probability; //recombination probability (#if ONECHROM == 0)

	int nL; //number of neutral loci
	//Deleterious mutations
	int loadEffect; //0 = affects offsping survival; 1 = affects fecundity and fertilization or mating probability in males
	int initial_nMut; //initial no. of deleterious mutations
	double R; //genome map length (see Roze & Rousset 2009, JEB) - corresponds to recombination rate
	double Ud; //mutation rate for mildly deleterious mutations / diploid genome / generation (see Spigler et al. 2016, Evolution)
	double Ul; //mutation rate for lethal mutations / diploid genome / generation (see Spigler et al. 2016, Evolution)
	double Ub; //mutation rate for back mutations / diploid genome / generation 
	double mean_sd; //mean selection coefficient for mildly deleterious mutations
	double mean_hd; //mean dominance coefficient for mildly deleterious mutations
	double sl; //selection coefficient for lethal mutations (fix - not sampled from a distribution)	
	double hl; //dominance coefficient for lethal mutations
	double sb; //selection coefficient for beneficial mutations

	//Reproduction & survival 
	bool c_f; //cost on female multiple mating
	double fec; //mean fecundity
	double omega_f; //(omega^2) strenght of selection on female multiple mating
	
	//Dispersal
	bool dispersal;
	bool postMating; // 0 = pre-mating dispersal; 1 = post-mating dispersal
	double dispCost; //paid in terms of survival

	//Mating system
	int matSys; //type of mating system: 0 = random polyginandorus; 1 = full polygyny (1 male for all females); 
				//2 = full promiscuity (1 random male per egg); 3 = monogamy

	void outPara(string dir); //Parameter output
};

