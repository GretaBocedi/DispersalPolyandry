#include "Parameters.h"

Parameters::Parameters() {
	//Simulation
	SimNr = 5;
	rep = 10; //replicates
	gen = 20001; //generations
	out_int = 50; //generation interval for output
	out_start = 0; //output start generation
	ind_interval = 500; //generations interval for individual outputs 
	indRand_interval = 500;
	PopMut_interval = 2500;
	out_mutations = 0;

	//Landscape
	spat_heterog = false; //spatial heterogeneity in K
	resol = 100; //resolution (m)
	x_max = 10;
	y_max = 10;
	cutoff = 40; //distance from the front after which rows should be cut off
	prop_suitable = 1.0; //proportion of suitable cells 
	K = 50.0; //carrying capacity (individuals/ha)
	Kmin = 2.0; //minimum K (for spatial heterogeneity)
	Kmax = 50.0; //max K

	exp_start = 50000; //start of range expansion
	start_dispEvol = 0; //start traits evolution
	max_age = 5; //max. age of a population to be considered part of the front

	//temporal environmental stochasticity
	env_stoch = false; 
	estd = 0.2; //amplitude for time series in K (Kmin will be zero; K needs to be set to mean K)
	ac = 0.0; //temporal autocorrelation for time series (Kmin will be zero; K needs to be set to mean K)
	p_cell_ext = 0.0; //probaility of local extinction

	//Initialisation
	min_seedX = 0;
	min_seedY = 0;
	max_seedX = x_max;
	max_seedY = y_max;

	//Traits
	sexDisp = true; //sex-biased dispersal
	nearest = false; //nearest-neighbour dispersal
	L = 1; //number of loci for each trait
	dispEvol = 1; //0 = no dispersal evolution; 1 = emigration p. evolving; 2 = disp. distance evolving 3 = both evolving
	polyEvol = 0; //0 = no polyandry evolution; 1 = polyandry evolving
	genot_mean = new double[5] {3.0, 0.05, 0.05, 200.0, 200.0}; // initial genotypic mean for each trait (a, ep_f, ep_m, dist_f, dist_m)
	genot_std = new double[5]{ 0.5, 0.1, 0.1, 0.00000000001, 0.00000000001 }; // initial genotypic standard deviation for each trait
	genot_std2 = new double[5]{ 0.1, 0.1, 0.1, 0.00000000001, 0.00000000001 }; //Needed to set mu_std
	mu_std = new double[5]{ 0.00000000001, 0.00000000001,0.00000000001, 0.00000000001, 0.00000000001 }; //DO NOT CHANGE
	mu = 0.001; // haploid per locus mutation rate for each trait
	mu_neutral = 0.001; // haploid per locus mutation probability for neutral loci

	rec_probability = 0.5; //recombination probability for adaptive loci (#if ONECHROM == 0)
	
	nL = 500; //number of neutral loci

	//Deleterious mutations
	loadEffect = 1; //0 = affects offsping survival; 1 = affects fecundity and fertilization or mating probability in males; -9 = fix inbreeding depression hard-coded in PolyDisp.cpp
	initial_nMut = 0;
	R = 10.0; //genome map length (see Roze & Rousset 2009, JEB) - corresponds to recombination rate
	Ud = 1.0; //mutation rate for mildly deleterious mutations / diploid genome / generation (see Spigler et al. 2016, Evolution)
	Ul = 0.2; //mutation rate for lethal mutations / diploid genome / generation (see Spigler et al. 2016, Evolution)
	Ub = 0.00000000001; //mutation rate for back mutations
	mean_sd = 0.05; //mean selection coefficient for mildly deleterious mutations
	mean_hd = 0.3; //mean dominance coefficient for mildly deleterious mutations
	sl = 1.0; //selection coefficient for lethal mutations (fix - not sampled from a distribution)
	hl = 0.02; //mean dominance coefficient for lethal mutations
	sb = -0.005; 

	//Reproduction & survival 
	c_f = 0; //cost on female multiple mating
	fec = 12.0; //mean fecundity
	omega_f = 700; //(omega^2) strenght of selection on female multiple mating

	//Dispersal
	dispersal = true;
	postMating = false; // 0 = pre-mating dispersal; 1 = post-mating dispersal
	dispCost = 0.6; //paid in terms of survival

	//Mating system
	matSys = 0; //type of mating system: 0 = random polyginandorus; 1 = full polygyny (1 male for all females); 
				//2 = full promiscuity (1 random male per egg); 3 = monogamy


}
//------------------------------
void Parameters::outPara(string name) {

	ofstream out;
	
	out.open(name.c_str());

	//Simulation
	out << "SimNr\t" << SimNr << endl;
	out << "rep\t" << rep << endl; //replicates
	out << "gen\t" << gen << endl; //generations
	out << "out_int\t" << out_int << endl;
	out << "out_start\t" << out_start << endl;
	out << "out_ind_interval\t" <<ind_interval << endl;
	out << "indRand_interval\t" << indRand_interval << endl;
	out << "PopMut_interval\t" << PopMut_interval << endl;
	out << "out_mutations\t" << out_mutations << endl;

	//Landscape
	out << "envir_stochasticity\t" << env_stoch << endl;
	out << "spatial_heterogenity_in_K\t" << spat_heterog << endl;
	out << "resolution\t" << resol << endl; //resolution (m)
	out << "x_max\t" << x_max << endl;
	out << "y_max\t" << y_max << endl;
	out << "cut_off\t" << cutoff << endl;
	out << "proportion_suitable_cells\t" << prop_suitable << endl; //proportion of suitable cells 	
	out << "K\t" << K << endl; //carrying capacity (individuals/ha)
	out << "Kmin\t" << Kmin << endl; //min carrying capacity for spatial heterogeneity
	out << "Kmax\t" << Kmax << endl; //max carrying capacity 
	out << "p_cell_ext\t" << p_cell_ext << endl;
	out << "expansion_start\t" << exp_start << endl;
	out << "dispersal_starts_evolving_at\t" << start_dispEvol << endl;
	out << "max_age\t" << max_age << endl;
	out << "temp_stochasticity_amplitude\t" << estd << endl;
	out << "temp_autocorrelation\t" << ac << endl;

	//Initialisation
	out << "min_seedX\t" << min_seedX << endl;
	out << "min_seedY\t" << min_seedY << endl;
	out << "max_seedX\t" << max_seedX << endl;
	out << "max_seedY\t" << max_seedY << endl;

	//Traits
	out << "sex-dependent_dispersal\t" << sexDisp << endl;
	out << "nearest-neigh\t" << nearest << endl;
	out << "nLoci_per_trait\t" << L << endl;
	out << "dispersal_evolving\t" << dispEvol << endl;
	out << "polyandry_evolving\t" << polyEvol << endl;
	out << "mean_poly\t" << genot_mean[0] << endl;
	out << "mean_emig_p_f\t" << genot_mean[1] << endl;
	out << "mean_emig_p_m\t" << genot_mean[2] << endl;
	out << "mean_disp_dist_f\t" << genot_mean[3] << endl;
	out << "mean_disp_dist_m\t" << genot_mean[4] << endl;
	out << "sd_poly\t" << genot_std[0] << endl;
	out << "sd_emig_p_f\t" << genot_std[1] << endl;
	out << "sd_emig_p_m\t" << genot_std[2] << endl;
	out << "sd_disp_dist_f\t" << genot_std[3] << endl;
	out << "sd_disp_dist_m\t" << genot_std[4] << endl;
	out << "mu\t" << mu << endl;
	out << "mu_neutral\t" << mu_neutral << endl;
	out << "mut_sd_poly\t" << mu_std[0] << endl;
	out << "mut_sd_emig_p_f\t" << mu_std[1] << endl;
	out << "mut_sd_emig_p_m\t" << mu_std[2] << endl;
	out << "mut_sd_disp_dist_f\t" << mu_std[3] << endl;
	out << "mut_sd_disp_dist_m\t" << mu_std[4] << endl;
	out << "recombination\t" << rec_probability << endl;

	out << "nNeutral_loci\t" << nL << endl;
	//Deleterious mutations
	out << "load_effect\t" << loadEffect << endl;
	out << "initial_Nmut\t" << initial_nMut << endl;
	out << "genome_map_length\t" << R << endl;
	out << "Ud\t" << Ud << endl;
	out << "Ul\t" << Ul << endl;
	out << "Ub\t" << Ub << endl;
	out << "mean_selection_mild_mut\t" << mean_sd << endl;
	out << "mean_dominance_mild_mut\t" << mean_hd << endl;
	out << "selection_lethal_mut\t" << sl << endl;
	out << "dominance_lethal_mut\t" << hl << endl;
	out << "selection_beneficial\t" << sb << endl;

	//Reproduction & survival 
	out << "fecundity\t" << fec << endl;
	out << "cost_polyandry\t" << c_f << endl;
	out << "omega_squared_polyandry\t" << omega_f << endl;
	out << "cost_dispersal\t" << dispCost << endl;
	out << "mating_system\t" << matSys << endl;

	out.close();
}
//------------------------------
Parameters::~Parameters()
{
}
