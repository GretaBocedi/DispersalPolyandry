#include "Individuals.h"

std::random_device rd2;
std::mt19937 gen(rd2());

std::bernoulli_distribution findhom(0.5); //for sampling the homologue


Individuals::Individuals(Parameters para, bool sx, int xx, int yy) {
	dispersed = false;
	alive = true;
	reproduce = true;
	sex = sx;
	x = xx; y = yy;
	xnatal = xx; ynatal = yy;
	noffs = 0;
	n_mates = 0;
	w = 1.0;
	h = 0;
	traits.initialise(para);
}
//--------------------------------------------------------------------------
Individuals::~Individuals() {
	
}
//--------------------------------------------------------------------------
void Individuals::deleteInd(void) {
	traits.deleteTraits();
#if MUTATIONLOAD
	chromo.deleteChromo();
#endif
	if (!mates.empty()) mates.clear();
	if (!markers.empty()) markers.clear();
}
//--------------------------------------------------------------------------

void Individuals::initialise(double kk, Parameters para, std::normal_distribution<> NormEm,  std::normal_distribution<> NormDist,
	std::normal_distribution<> NormPoly, std::gamma_distribution<> finds, std::uniform_real_distribution<> neutr,
	std::uniform_real_distribution<> findpos) {

	int hom;
	double pos;
	double ss;
	double hh;

	//initialise neutral markers
	for (int i = 0; i < (para.nL * 2); i++) markers.push_back(neutr(gen));


#if MUTATIONLOAD
	if (para.initial_nMut > 0) {

		for (int i = 0; i < para.initial_nMut; i++) {
			//sample homologue
			hom = findhom(gen);
			//sample selection coefficient
			ss = finds(gen);
			//sample dominance coefficient
			std::uniform_real_distribution<> findh(0.0, std::exp(-kk * ss));
			hh = findh(gen);

			//sample position
			pos = findpos(gen);
			//add mutation
			w *= chromo.addDelMutation(hom, pos, ss, hh);
		}
	}
#endif
	
	if (para.polyEvol) {
		for (int i = 0; i < (2 * para.L); i++) {
			traits.a_mat[i] = NormPoly(gen);
			*traits.g_a_mat += traits.a_mat[i];
		}
		*traits.p_a_mat = *traits.g_a_mat;
		if (*traits.p_a_mat < 0.0) *traits.p_a_mat = 0.0;
	}

	if (para.dispEvol > 0) {
		//if dispersal is not sex limited use only the female trait
		switch (para.dispEvol) {
		case 1:
			for (int i = 0; i < (2 * para.L); i++) {
				traits.ep_f[i] = NormEm(gen);
				*traits.g_ep_f += traits.ep_f[i];
			}
			*traits.p_ep_f = *traits.g_ep_f;
			if (*traits.p_ep_f < 0.0) *traits.p_ep_f = 0.0;
			if (*traits.p_ep_f > 1.0) *traits.p_ep_f = 1.0;
			break;
		case 2:
			for (int i = 0; i < (2 * para.L); i++) {
				traits.dist_f[i] = NormDist(gen);
				*traits.g_dist_f += traits.dist_f[i];
			}
			*traits.p_dist_f = *traits.g_dist_f;
			if (*traits.p_dist_f < 0.0) *traits.p_dist_f = 0.0;
			break;
		case 3:
			for (int i = 0; i < (2 * para.L); i++) {
				traits.ep_f[i] = NormEm(gen);
				*traits.g_ep_f += traits.ep_f[i];
				traits.dist_f[i] = NormDist(gen);
				*traits.g_dist_f += traits.dist_f[i];
			}
			*traits.p_ep_f = *traits.g_ep_f;
			if (*traits.p_ep_f < 0.0) *traits.p_ep_f = 0.0;
			if (*traits.p_ep_f > 1.0) *traits.p_ep_f = 1.0;
			*traits.p_dist_f = *traits.g_dist_f;
			if (*traits.p_dist_f < 0.0) *traits.p_dist_f = 0.0;
			break;
		}
	}
}
//---------------------------------------------------------------------------------------------
//Initialised when dispersal is sex-limited
void Individuals::initialise(double kk, Parameters para, std::normal_distribution<> NormEm, std::normal_distribution<> NormDist,
	std::normal_distribution<> NormPoly, std::gamma_distribution<> finds, std::uniform_real_distribution<> neutr,
	std::uniform_real_distribution<> findpos, std::normal_distribution<> NormEmM, std::normal_distribution<> NormDistM) {

	int hom;
	double pos;
	double ss;
	double hh;

	//initialise neutral markers
	for (int i = 0; i < (para.nL * 2); i++) markers.push_back(neutr(gen));


#if MUTATIONLOAD
	if (para.initial_nMut > 0) {

		for (int i = 0; i < para.initial_nMut; i++) {
			//sample homologue
			hom = findhom(gen);
			//sample selection coefficient
			ss = finds(gen);
			//sample dominance coefficient
			std::uniform_real_distribution<> findh(0.0, std::exp(-kk * ss));
			hh = findh(gen);

			//sample position
			pos = findpos(gen);
			//add mutation
			w *= chromo.addDelMutation(hom, pos, ss, hh);
		}
	}
#endif

	if (para.polyEvol) {
		for (int i = 0; i < (2 * para.L); i++) {
			traits.a_mat[i] = NormPoly(gen);
			*traits.g_a_mat += traits.a_mat[i];
		}
		*traits.p_a_mat = *traits.g_a_mat;
		if (*traits.p_a_mat < 0.0) *traits.p_a_mat = 0.0;
	}

	if (para.dispEvol > 0) {
		switch (para.dispEvol) {
		case 1:
			//female emigration probability
			for (int i = 0; i < (2 * para.L); i++) {
				traits.ep_f[i] = NormEm(gen);
				*traits.g_ep_f += traits.ep_f[i];
			}
			*traits.p_ep_f = *traits.g_ep_f;
			if (*traits.p_ep_f < 0.0) *traits.p_ep_f = 0.0;
			if (*traits.p_ep_f > 1.0) *traits.p_ep_f = 1.0;
			//male emigration probability
			for (int i = 0; i < (2 * para.L); i++) {
				traits.ep_m[i] = NormEmM(gen);
				*traits.g_ep_m += traits.ep_m[i];
			}
			*traits.p_ep_m = *traits.g_ep_m;
			if (*traits.p_ep_m < 0.0) *traits.p_ep_m = 0.0;
			if (*traits.p_ep_m > 1.0) *traits.p_ep_m = 1.0;
			break;
		case 2:
			//female dispersal distance
			for (int i = 0; i < (2 * para.L); i++) {
				traits.dist_f[i] = NormDist(gen);
				*traits.g_dist_f += traits.dist_f[i];
			}
			*traits.p_dist_f = *traits.g_dist_f;
			if (*traits.p_dist_f < 0.0) *traits.p_dist_f = 0.0;
			//male dispersal distance
			for (int i = 0; i < (2 * para.L); i++) {
				traits.dist_m[i] = NormDistM(gen);
				*traits.g_dist_m += traits.dist_m[i];
			}
			*traits.p_dist_m = *traits.g_dist_m;
			if (*traits.p_dist_m < 0.0) *traits.p_dist_m = 0.0;
			break;
		case 3:
			//female emigration probablilty and dispersal distance
			for (int i = 0; i < (2 * para.L); i++) {
				traits.ep_f[i] = NormEm(gen);
				*traits.g_ep_f += traits.ep_f[i];
				traits.dist_f[i] = NormDist(gen);
				*traits.g_dist_f += traits.dist_f[i];
			}
			*traits.p_ep_f = *traits.g_ep_f;
			if (*traits.p_ep_f < 0.0) *traits.p_ep_f = 0.0;
			if (*traits.p_ep_f > 1.0) *traits.p_ep_f = 1.0;
			*traits.p_dist_f = *traits.g_dist_f;
			if (*traits.p_dist_f < 0.0) *traits.p_dist_f = 0.0;
			//male emigration probaibilty and dispersal distance
			for (int i = 0; i < (2 * para.L); i++) {
				traits.ep_m[i] = NormEmM(gen);
				*traits.g_ep_m += traits.ep_m[i];
				traits.dist_m[i] = NormDistM(gen);
				*traits.g_dist_m += traits.dist_m[i];
			}
			*traits.p_ep_m = *traits.g_ep_m;
			if (*traits.p_ep_m < 0.0) *traits.p_ep_m = 0.0;
			if (*traits.p_ep_m > 1.0) *traits.p_ep_m = 1.0;
			*traits.p_dist_m = *traits.g_dist_m;
			if (*traits.p_dist_m < 0.0) *traits.p_dist_m = 0.0;
			break;
		}
	}
}
//--------------------------------------------------------------------------
void Individuals::neutral_mutation(int nmut, std::uniform_int_distribution<> findpos, std::uniform_real_distribution<> mutate) {
	int allele;

	for (int i = 0; i < nmut; i++) {
		allele = findpos(gen);
		//if locus was homozygote, discount for that
		if (allele % 2 == 0.0) { //1st allele
			if (markers[allele] == markers[allele + 1]) h--;
		}
		else { //2nd allele
			if (markers[allele - 1] == markers[allele]) h--;
		}
		//assume that you won't get an homozygote from continuous mutation 
		markers[allele] = mutate(gen);
	}
}
//--------------------------------------------------------------------------
#if MUTATIONLOAD
void Individuals::benef_mutation(int nmut, double s, std::uniform_real_distribution<> findpos) {
	int hom;
	double pos;
	double ss;
	std::bernoulli_distribution findhom(0.5); //for sampling the homologue

	for (int i = 0; i < nmut; i++) {
		//sample homologue
		hom = findhom(gen);
		//sample selection coefficient
		ss = s;
		//sample position
		pos = findpos(gen);
		//add mutation
		w *= chromo.addDelMutation(hom, pos, ss, 1.0); //assume complete dominance of beneficial mutations
	}
}
//--------------------------------------------------------------------------
void Individuals::back_mutation(int nmut)
{
	int rdn;
	std::map<double, mutation, std::less<double>>::iterator it;

	std::bernoulli_distribution hom(0.5);

	for (int i = 0; i < nmut; i++) {
		std::uniform_int_distribution<> samp(0, chromo.mutations.size() - 1);
		rdn = samp(gen);
		it = chromo.mutations.begin();
		std::advance(it, rdn);

		if (it->second.homol == 2) { //mutation is homozygote
			it->second.homol = hom(gen);
			chromo.Nho--;		
			//update fitness
			w /= (1.0 - it->second.s);
			w *= (1.0 - it->second.h * it->second.s);
		}
		else { //mutation is heterozygote
			chromo.mutations.erase(it);
			chromo.nMut--;
			//update fitness
			w /= (1.0 - it->second.h * it->second.s);
		}
	}
}
//--------------------------------------------------------------------------
void Individuals::delet_mutation(int nmut, double kk, std::uniform_real_distribution<> findpos, 
	 std::gamma_distribution<> finds) {	
	int hom;
	double pos;
	double ss;
	double hh;


	for (int i = 0; i < nmut; i++) {
		//sample homologue
		hom = findhom(gen);
		//sample selection coefficient
		ss = finds(gen);
		//sample dominance coefficient
		std::uniform_real_distribution<> findh(0.0, std::exp(-kk * ss));
		hh = findh(gen);
		//sample position
		pos = findpos(gen);
		//add mutation
		w *= chromo.addDelMutation(hom, pos, ss, hh);
	}
}
//--------------------------------------------------------------------------
void Individuals::delet_mutation(int nmut, std::uniform_real_distribution<> findpos, 
	double ss, double hh) {
	int hom;
	double pos;

	for (int i = 0; i < nmut; i++) {
		hom = findhom(gen);
		pos = findpos(gen);
		w *= chromo.addDelMutation(hom, pos, ss, hh);
	}
}
#endif
//--------------------------------------------------------------------------
void Individuals::traits_mutation(int trait, int allele,  std::normal_distribution<> Norm) {

	switch (trait) {
	case 0: //re-mating slope (polyandry)
		*traits.g_a_mat -= traits.a_mat[allele];
		traits.a_mat[allele] += Norm(gen);
		*traits.g_a_mat += traits.a_mat[allele];
		*traits.p_a_mat = *traits.g_a_mat;
		if (*traits.p_a_mat < 0.0) *traits.p_a_mat = 0.0;
		break;
	case 1: //female emigration probability
		*traits.g_ep_f -= traits.ep_f[allele];
		traits.ep_f[allele] += Norm(gen);
		*traits.g_ep_f += traits.ep_f[allele];
		*traits.p_ep_f = *traits.g_ep_f;
		if (*traits.p_ep_f < 0.0) *traits.p_ep_f = 0.0;
		if (*traits.p_ep_f > 1.0) *traits.p_ep_f = 1.0;
		break;
	case 2: //female dispersal distance
		*traits.g_dist_f -= traits.dist_f[allele];
		traits.dist_f[allele] += Norm(gen);
		*traits.g_dist_f += traits.dist_f[allele];
		*traits.p_dist_f = *traits.g_dist_f;
		if (*traits.p_dist_f < 0.0) *traits.p_dist_f = 0.0;
		break;
	case 3: //male emigration probability
		*traits.g_ep_m -= traits.ep_m[allele];
		traits.ep_m[allele] += Norm(gen);
		*traits.g_ep_m += traits.ep_m[allele];
		*traits.p_ep_m = *traits.g_ep_m;
		if (*traits.p_ep_m < 0.0) *traits.p_ep_m = 0.0;
		if (*traits.p_ep_m > 1.0) *traits.p_ep_m = 1.0;
		break;
	case 4: //male dispersal distance
		*traits.g_dist_m -= traits.dist_m[allele];
		traits.dist_m[allele] += Norm(gen);
		*traits.g_dist_m += traits.dist_m[allele];
		*traits.p_dist_m = *traits.g_dist_m;
		if (*traits.p_dist_m < 0.0) *traits.p_dist_m = 0.0;
		break;
	}
}
//--------------------------------------------------------------------------
void Individuals::outMut(int g, std::ofstream *out) {

	map<double, mutation>::iterator iter;

#if MUTATIONLOAD
	for (iter = chromo.mutations.begin(); iter != chromo.mutations.end(); iter++) {
		*out << g << "\t" << iter->second.s << "\t" << iter->second.h << "\t";
		if (iter->second.homol < 2) *out << 0;
		else *out << 1;
		*out << endl;
	}
#endif

}
//--------------------------------------------------------------------------
void Individuals::outRanInd(int outb, int r, int g, std::ofstream * out)
{
#if MUTATIONLOAD
	*out << r << "\t" << g << "\t" << x << "\t" << y << "\t" << outb << "\t" << w << "\t" << h << "\t" << chromo.nMut << "\t" << chromo.Nho << endl;
#endif
}
//--------------------------------------------------------------------------
void Individuals::outInd(Parameters para, int r, int g, std::ofstream *out) {
	*out << r << "\t" << g << "\t" << x << "\t" << y << "\t" << sex << "\t";
	*out << alive << "\t" << reproduce << "\t" << dispersed;
	*out << "\t" << n_mates << "\t" << w << "\t" << h;
#if MUTATIONLOAD
	*out << "\t" << chromo.nMut << "\t" << chromo.Nho;
#endif
	if (para.polyEvol) *out << "\t" << *traits.g_a_mat << "\t" << *traits.p_a_mat;
	if (para.dispEvol > 0) {
		switch (para.dispEvol) {
		case 1:
			if (para.sexDisp) {
				if (sex) *out << "\t" << *traits.g_ep_f << "\t" << *traits.p_ep_f;
				else *out << "\t" << *traits.g_ep_m << "\t" << *traits.p_ep_m;
			}
			else *out << "\t" << *traits.g_ep_f << "\t" << *traits.p_ep_f;
			break;
		case 2:
			if (para.sexDisp) {
				if (sex) *out << "\t" << *traits.g_dist_f << "\t" << *traits.p_dist_f;
				else *out << "\t" << *traits.g_dist_m << "\t" << *traits.p_dist_m;
			}
			else *out << "\t" << *traits.g_dist_f << "\t" << *traits.p_dist_f;
			break;
		case 3:
			if (para.sexDisp) {
				if (sex) *out << "\t" << *traits.g_ep_f << "\t" << *traits.p_ep_f << "\t" << *traits.g_dist_f << "\t" << *traits.p_dist_f;
				else *out << "\t" << *traits.g_ep_m << "\t" << *traits.p_ep_m << "\t" << *traits.g_dist_m << "\t" << *traits.p_dist_m;
			}
			else *out << "\t" << *traits.g_ep_f << "\t" << *traits.p_ep_f << "\t" << *traits.g_dist_f << "\t" << *traits.p_dist_f;
			break;
		}
	}
	*out << endl;
}

//--------------------------------------------------------------------------
void Individuals::outOffs(Parameters para, int r, int g, std::ofstream* out) {
	*out << r << "\t" << g << "\t" << x << "\t" << y << "\t" << sex << "\t";
	*out << alive << "\t" << w << "\t" << h;
#if MUTATIONLOAD
	*out << "\t" << chromo.nMut << "\t" << chromo.Nho;
#endif
	*out << endl;
}

