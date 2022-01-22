#pragma once

#include <map>
#include <set>

#include <stdio.h>
#include <fstream>
#include <iostream>

#include "Parameters.h"


struct mutation {

	int homol; //0 = only on homologue 1; 1 = only on homologue 2; 2 = on both homologues (hence homozygote)
	double s; //selection coefficient
	double h; //dominance coefficient
};

//map containing all mutations: map<position, mutation> 
typedef std::map<double, mutation, std::less<double>> MapMuts;

typedef std::map<double, int, std::less<double>> MapPos;

//map containing continuous allele coding for a trait: map<position, allele>
typedef std::map<double, double, std::less<double>> MapTrait;


class Chromosome {
public:
	Chromosome();
	~Chromosome();
	int nMut; //no. of deleterious mutations
	int Nho; //number of homozygote mutations

	MapMuts mutations; //list of all deleterious mutations, their coefficients and whether they are homo or heterozygote 
	
	double addDelMutation(
		int, //homologue
		double, //position
		double, //s
		double //h
		);
	void deleteChromo();

#if ONECHROM	
	MapTrait a_mat1, a_mat2; //Maps of allelic values on the two homologues for re-mating slope (probability of mating = exp(-a_mat * n_mates))
	MapTrait ep_f1, ep_f2; //Maps of allelic values on the two homologues for female emigration probability
	MapTrait ep_m1, ep_m2; //Maps of allelic values on the two homologues for male emigration probability
	MapTrait dist_f1, dist_f2; //Maps of allelic values on the two homologues for female mean dispersal distance
	MapTrait dist_m1, dist_m2; //Maps of allelic values on the two homologues for male  mean dispersal distance
#endif
private:
};