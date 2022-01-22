#include "Population.h"
#include "Individuals.h"


Population::Population(int xx, int yy) {
	age = 0;
	x = xx;
	y = yy;
	N = 0;
	Noffs = 0;
	Nf = 0;
	Nm = 0;
	Foffs = 0;
	Moffs = 0;
	DelMut = 0;
	meanDelMut = 0.0;
	DelMutSd = 0.0; 
	W = 0.0;
	Wsd = 0.0;
	Homoz = 0.0; Hsd = 0.0;
	mMates = 0.0; sMates = 0.0;
	mA = 0.0; mFep = 0.0; mMep = 0.0; mFdist = 0.0; mMdist = 0.0; 
	sA = 0.0; sFep = 0.0; sMep = 0.0; sFdist = 0.0; sMdist = 0.0; 
	mgA = 0.0; mgFep = 0.0; mgMep = 0.0; mgFdist = 0.0; mgMdist = 0.0;
	sgA = 0.0; sgFep = 0.0; sgMep = 0.0; sgFdist = 0.0; sgMdist = 0.0;

	if (!females.empty()) females.clear();
	if (!males.empty()) males.clear();
	if (!Jfemales.empty()) Jfemales.clear();
	if (!Jmales.empty()) Jmales.clear();
	if (!tmp_females.empty()) tmp_females.clear();
	if (!tmp_males.empty()) tmp_males.clear();
	if (!popMuts.empty()) popMuts.clear();
}
//----------------------------------------
Population::~Population() {
	if (!females.empty()) females.clear();
	if (!males.empty()) males.clear();
	if (!Jfemales.empty()) Jfemales.clear();
	if (!Jmales.empty()) Jmales.clear();
	if (!tmp_females.empty()) tmp_females.clear();
	if (!tmp_males.empty()) tmp_males.clear();
}
//----------------------------------------
void Population::initialise_pop(double k, double kk, Parameters para, std::normal_distribution<> NormEm, std::normal_distribution<> NormDist,
	std::normal_distribution<> NormPoly, std::gamma_distribution<> finds, std::uniform_real_distribution<> neutr, 
	std::uniform_real_distribution<> findpos) {

	N = (int)k;

	//initialise females
	for (int i = 0; i < (int)(N/2); i++) {
		females.push_back(Individuals(para, true, x, y));	
		females[i].initialise(kk, para, NormEm, NormDist, NormPoly, finds, neutr, findpos);
		Nf++;
		W += 1.0; //at initialisation every individual has viability = 1
	}
	//initialise males
	for (int i = 0; i < (int)N/2; i++) {
		males.push_back(Individuals(para, false, x, y));
		males[i].initialise(kk, para, NormEm, NormDist, NormPoly, finds, neutr, findpos);
		Nm++;
		W += 1.0; //at initialisation every individual has viability = 1
	}
}
//----------------------------------------
// Initialise pop when dispersal is sex-limited
void Population::initialise_pop(double k, double kk, Parameters para, std::normal_distribution<> NormEm, std::normal_distribution<> NormDist,
	std::normal_distribution<> NormPoly, std::gamma_distribution<> finds, std::uniform_real_distribution<> neutr,
	std::uniform_real_distribution<> findpos, std::normal_distribution<> NormEmM, std::normal_distribution<> NormDistM) {

	N = (int)k;

	//initialise females
	for (int i = 0; i < (int)(N / 2); i++) {
		females.push_back(Individuals(para, true, x, y));
		females[i].initialise(kk, para, NormEm, NormDist, NormPoly, finds, neutr, findpos, NormEmM, NormDistM);
		Nf++;
		W += 1.0; //at initialisation every individual has viability = 1
	}
	//initialise males
	for (int i = 0; i < (int)N / 2; i++) {
		males.push_back(Individuals(para, false, x, y));
		males[i].initialise(kk, para, NormEm, NormDist, NormPoly, finds, neutr, findpos, NormEmM, NormDistM);
		Nm++;
		W += 1.0; //at initialisation every individual has viability = 1
	}
}
//----------------------------------------
void Population::computeSums(Parameters para, Individuals ind) {
#if MUTATIONLOAD
	DelMut += ind.chromo.nMut;
	DelMutSd += (double)ind.chromo.nMut * (double)ind.chromo.nMut; //sum of squares
	W += ind.w;
	Wsd += ind.w * ind.w; //sum of squares
#endif
	//neutral homozygosity
	Homoz += (double)ind.h / (double)para.nL;
	Hsd += ((double)ind.h / (double)para.nL) * ((double)ind.h / (double)para.nL);//sum of squares

	if (ind.sex) {
		mMates += ind.n_mates;
		sMates += (ind.n_mates * ind.n_mates);
	}
	if (para.polyEvol && ind.sex) {
		mA += *ind.traits.p_a_mat;
		sA += (*ind.traits.p_a_mat) * (*ind.traits.p_a_mat);
		mgA += *ind.traits.g_a_mat;
		sgA += (*ind.traits.g_a_mat) * (*ind.traits.g_a_mat);
	}
	if (para.dispEvol > 0) {
		if (para.sexDisp) {
			if (ind.sex) {
				switch (para.dispEvol) {
				case 1:
					mFep += *ind.traits.p_ep_f; sFep += (*ind.traits.p_ep_f) * (*ind.traits.p_ep_f);
					mgFep += *ind.traits.g_ep_f; sgFep += (*ind.traits.g_ep_f) * (*ind.traits.g_ep_f);
					break;
				case 2:
					mFdist += *ind.traits.p_dist_f; sFdist += (*ind.traits.p_dist_f) * (*ind.traits.p_dist_f);
					mgFdist += *ind.traits.g_dist_f; sgFdist += (*ind.traits.g_dist_f) * (*ind.traits.g_dist_f);
					break;
				case 3:
					mFep += *ind.traits.p_ep_f; sFep += (*ind.traits.p_ep_f) * (*ind.traits.p_ep_f);
					mFdist += *ind.traits.p_dist_f; sFdist += (*ind.traits.p_dist_f) * (*ind.traits.p_dist_f);
					mgFep += *ind.traits.g_ep_f; sgFep += (*ind.traits.g_ep_f) * (*ind.traits.g_ep_f);
					mgFdist += *ind.traits.g_dist_f; sgFdist += (*ind.traits.g_dist_f) * (*ind.traits.g_dist_f);
					break;
				}
			}
			else {
				switch (para.dispEvol) {
				case 1:
					mMep += *ind.traits.p_ep_m; sMep += (*ind.traits.p_ep_m) * (*ind.traits.p_ep_m);
					mgMep += *ind.traits.g_ep_m; sgMep += (*ind.traits.g_ep_m) * (*ind.traits.g_ep_m);
					break;
				case 2:
					mMdist += *ind.traits.p_dist_m; sMdist += (*ind.traits.p_dist_m) * (*ind.traits.p_dist_m);
					mgMdist += *ind.traits.g_dist_m; sgMdist += (*ind.traits.g_dist_m) * (*ind.traits.g_dist_m);
					break;
				case 3:
					mMep += *ind.traits.p_ep_m; sMep += (*ind.traits.p_ep_m) * (*ind.traits.p_ep_m);
					mMdist += *ind.traits.p_dist_m; sMdist += (*ind.traits.p_dist_m) * (*ind.traits.p_dist_m);
					mgMep += *ind.traits.g_ep_m; sgMep += (*ind.traits.g_ep_m) * (*ind.traits.g_ep_m);
					mgMdist += *ind.traits.g_dist_m; sgMdist += (*ind.traits.g_dist_m) * (*ind.traits.g_dist_m);
					break;
				}
			}
		}
		else {
			switch (para.dispEvol) {
			case 1:
				mFep += *ind.traits.p_ep_f; sFep += (*ind.traits.p_ep_f) * (*ind.traits.p_ep_f);
				mgFep += *ind.traits.g_ep_f; sgFep += (*ind.traits.g_ep_f) * (*ind.traits.g_ep_f);
				break;
			case 2:
				mFdist += *ind.traits.p_dist_f; sFdist += (*ind.traits.p_dist_f) * (*ind.traits.p_dist_f);
				mgFdist += *ind.traits.g_dist_f; sgFdist += (*ind.traits.g_dist_f) * (*ind.traits.g_dist_f);
				break;
			case 3:
				mFep += *ind.traits.p_ep_f; sFep += (*ind.traits.p_ep_f) * (*ind.traits.p_ep_f);
				mFdist += *ind.traits.p_dist_f; sFdist += (*ind.traits.p_dist_f) * (*ind.traits.p_dist_f);
				mgFep += *ind.traits.g_ep_f; sgFep += (*ind.traits.g_ep_f) * (*ind.traits.g_ep_f);
				mgFdist += *ind.traits.g_dist_f; sgFdist += (*ind.traits.g_dist_f) * (*ind.traits.g_dist_f);
				break;
			}
		}
	}
}
//----------------------------------------
void Population::computeStats(Parameters para) {
#if MUTATIONLOAD
	DelMutSd = std::sqrt((DelMutSd - ((double)DelMut * (double)DelMut) / (double)N) / (double)N);
	meanDelMut = (double)DelMut / (double)N;

	Wsd = std::sqrt((Wsd - (W * W) / (double)N) / (double)N);
	W /= (double)N;

	Hsd = std::sqrt((Hsd - (Homoz * Homoz) / (double)N) / (double)N);
	Homoz /= (double)N;
#endif
	if (Nf > 0) {
		sMates = std::sqrt((sMates - (mMates * mMates) / (double)Nf) / (double)Nf);
		mMates /= (double)Nf;
	}

	if (para.polyEvol && Nf > 0) {
		sA = std::sqrt((sA - (mA * mA) / (double)Nf) / (double)Nf);
		mA /= (double)Nf;

		sgA = std::sqrt((sgA - (mgA * mgA) / (double)Nf) / (double)Nf);
		mgA /= (double)Nf;
	}
	if (para.dispEvol > 0) {
		if (para.sexDisp) {
			switch (para.dispEvol) {
			case 1:
				if (Nf > 0) {
					sFep = std::sqrt((sFep - (mFep * mFep) / (double)Nf) / (double)Nf);
					mFep /= (double)Nf;
					sgFep = std::sqrt((sgFep - (mgFep * mgFep) / (double)Nf) / (double)Nf);
					mgFep /= (double)Nf;

				}
				if (Nm > 0) {
					sMep = std::sqrt((sMep - (mMep * mMep) / (double)Nm) / (double)Nm);
					mMep /= (double)Nm;
					sgMep = std::sqrt((sgMep - (mgMep * mgMep) / (double)Nm) / (double)Nm);
					mgMep /= (double)Nm;
				}
				break;
			case 2:
				if (Nf > 0) {
					sFdist = std::sqrt((sFdist - (mFdist * mFdist) / (double)Nf) / (double)Nf);
					mFdist /= (double)Nf;
					sgFdist = std::sqrt((sgFdist - (mgFdist * mgFdist) / (double)Nf) / (double)Nf);
					mgFdist /= (double)Nf;
				}
				if (Nm > 0) {
					sMdist = std::sqrt((sMdist - (mMdist * mMdist) / (double)Nm) / (double)Nm);
					mMdist /= (double)Nm;
					sgMdist = std::sqrt((sgMdist - (mgMdist * mgMdist) / (double)Nm) / (double)Nm);
					mgMdist /= (double)Nm;
				}
				break;
			case 3:
				if (Nf > 0) {
					sFep = std::sqrt((sFep - (mFep * mFep) / (double)Nf) / (double)Nf);
					mFep /= (double)Nf;
					sFdist = std::sqrt((sFdist - (mFdist * mFdist) / (double)Nf) / (double)Nf);
					mFdist /= (double)Nf;

					sgFep = std::sqrt((sgFep - (mgFep * mgFep) / (double)Nf) / (double)Nf);
					mgFep /= (double)Nf;
					sgFdist = std::sqrt((sgFdist - (mgFdist * mgFdist) / (double)Nf) / (double)Nf);
					mgFdist /= (double)Nf;
				}
				if (Nm > 0) {
					sMep = std::sqrt((sMep - (mMep * mMep) / (double)Nm) / (double)Nm);
					mMep /= (double)Nm;
					sMdist = std::sqrt((sMdist - (mMdist * mMdist) / (double)Nm) / (double)Nm);
					mMdist /= (double)Nm;

					sgMep = std::sqrt((sgMep - (mgMep * mgMep) / (double)Nm) / (double)Nm);
					mgMep /= (double)Nm;
					sgMdist = std::sqrt((sgMdist - (mgMdist * mgMdist) / (double)Nm) / (double)Nm);
					mgMdist /= (double)Nm;
				}
				break;
			}

		}
		else {
			switch (para.dispEvol) {
			case 1:
				sFep = std::sqrt((sFep - (mFep * mFep) / (double)N) / (double)N);
				mFep /= (double)N;
				sgFep = std::sqrt((sgFep - (mgFep * mgFep) / (double)N) / (double)N);
				mgFep /= (double)N;
				break;
			case 2:
				sFdist = std::sqrt((sFdist - (mFdist * mFdist) / (double)N) / (double)N);
				mFdist /= (double)N;
				sgFdist = std::sqrt((sgFdist - (mgFdist * mgFdist) / (double)N) / (double)N);
				mgFdist /= (double)N;
				break;
			case 3:
				sFep = std::sqrt((sFep - (mFep * mFep) / (double)N) / (double)N);
				mFep /= (double)N;
				sFdist = std::sqrt((sFdist - (mFdist * mFdist) / (double)N) / (double)N);
				mFdist /= (double)N;

				sgFep = std::sqrt((sgFep - (mgFep * mgFep) / (double)N) / (double)N);
				mgFep /= (double)N;
				sgFdist = std::sqrt((sgFdist - (mgFdist * mgFdist) / (double)N) / (double)N);
				mgFdist /= (double)N;
				break;
			}
		}
	}	
}
//----------------------------------------
void Population::set2zero(void) {
#if MUTATIONLOAD
	DelMut = 0;
	DelMutSd = 0.0;
	W = 0.0;
	Wsd = 0.0;
	Homoz = 0.0; Hsd = 0.0;
#endif
	mMates = 0.0; sMates = 0.0;
	mA = 0.0; mFep = 0.0; mMep = 0.0; mFdist = 0.0; mMdist = 0.0; 
	sA = 0.0; sFep = 0.0; sMep = 0.0; sFdist = 0.0; sMdist = 0.0; 

	mgA = 0.0; mgFep = 0.0; mgMep = 0.0; mgFdist = 0.0; mgMdist = 0.0;
	sgA = 0.0; sgFep = 0.0; sgMep = 0.0; sgFdist = 0.0; sgMdist = 0.0;
}
//----------------------------------------
void Population::deleteAdults(void) {
	vector<Individuals>::iterator iter;

	if (!females.empty()) {
		for (iter = females.begin(); iter != females.end(); iter++) iter->deleteInd();
		females.clear();
	}
	if (!males.empty()) {
		for (iter = males.begin(); iter != males.end(); iter++) iter->deleteInd();
		males.clear();
	}

	if (!tmp_females.empty()) tmp_females.clear();
	if (!tmp_males.empty()) tmp_males.clear();

}
//----------------------------------------
void Population::outPop(int r, int g, std::ofstream *out) {
	*out << r << "\t" << g << "\t" << x << "\t" << y << "\t" << age << "\t" << Nf << "\t" << Nm;
#if MUTATIONLOAD
	*out << "\t" << W << "\t" << Wsd << "\t";
	*out << Homoz << "\t" << Hsd << "\t" << meanDelMut << "\t" << DelMutSd;
#endif
	*out << endl;
}
//----------------------------------------
void Population::outTrait(Parameters para, int r, int g, std::ofstream *out) {
	*out << r << "\t" << g << "\t" << x << "\t" << y;
	*out << "\t" << mMates << "\t" << sMates;
	if (para.polyEvol) *out << "\t" << mA << "\t" << sA << "\t" << mgA << "\t" << sgA;
	if (para.dispEvol > 0) {
		switch (para.dispEvol) {
		case 1:
			*out << "\t" << mFep << "\t" << sFep << "\t" << mgFep << "\t" << sgFep;
			break;
		case 2:
			*out << "\t" << mFdist << "\t" << sFdist << "\t" << mgFdist << "\t" << sgFdist;
			break;
		case 3:
			*out << "\t" << mFep << "\t" << sFep << "\t" << mgFep << "\t" << sgFep;
			*out << "\t" << mFdist << "\t" << sFdist << "\t" << mgFdist << "\t" << sgFdist;
			break;
		}
		if (para.sexDisp) {
			switch (para.dispEvol) {
			case 1:
				*out << "\t" << mMep << "\t" << sMep << "\t" << mgMep << "\t" << sgMep;
				break;
			case 2:
				*out << "\t" << mMdist << "\t" << sMdist << "\t" << mgMdist << "\t" << sgMdist;
				break;
			case 3:
				*out << "\t" << mMep << "\t" << sMep << "\t" << mgMep << "\t" << sgMep;
				*out << "\t" << mMdist << "\t" << sMdist << "\t" <<  mgMdist << "\t" << sgMdist;
				break;
			}
		}
	}
	*out << endl;
}

void Population::addMutation(double pos, double ss, double hh)
{
	pop_muts muts;

	muts.h = hh;
	muts.s = ss;
	
	if (popMuts.find(pos) == popMuts.end()) {
		muts.count = 1;
		popMuts[pos] = muts;
	}
	else popMuts[pos].count++;

}

void Population::outMutations(int n, int r, int g, std::ofstream *out)
{
	map<double, pop_muts>::iterator iter;

	for (iter = popMuts.begin(); iter != popMuts.end(); iter++) {
		*out << r << "\t" << g << "\t" << x << "\t" << y << "\t";
		*out << iter->second.s << "\t" << iter->second.h << "\t";
		*out << (double)iter->second.count / (2.0 * (double)n) << endl;
	}
}



