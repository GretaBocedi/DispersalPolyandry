#include "PolyDisp.h"

//---------------------------------------------------------------------------
// Main function
#if CLUSTER
int main(int argc, char* argv[])
{
	// Get the current directory.
	char* buffer = getcwd(NULL, 0);
	dir = buffer;
	free(buffer);
	dir = dir + "/"; //Current directory path
	dirOut = dir + "Outputs/"; //Outpus folder path


	para.SimNr = std::atoi(argv[1]);
	para.rep = std::atoi(argv[2]);
	//para.prop_suitable = std::atof(argv[3]);
	//para.dispEvol = std::atoi(argv[3]);
	//para.polyEvol = std::atoi(argv[4]);
	/*para.genot_mean[1] = std::atof(argv[4]);
	para.genot_mean[2] = std::atof(argv[5]);*/
	para.dispCost = std::atof(argv[3]);
	/*para.Ul = std::atof(argv[7]);
	para.costIncrement = std::atof(argv[8]);
	para.Ud = std::atof(argv[9]);
	para.loadEffect = std::atoi(argv[10]);
	para.p_cell_ext = std::atof(argv[11]);*/
	para.c_f = std::atoi(argv[4]);
	para.omega_f = std::atof(argv[5]);

	RunModel();

	cout << "Simulation completed" << endl;

	return 0;
}
#else
int _tmain(int argc, _TCHAR* argv[])
{

	// Get the current directory.
	char* buffer = _getcwd(NULL, 0);
	dir = buffer;
	free(buffer);
	dir = dir + "\\"; //Current directory path
	dirOut = dir + "Outputs\\"; //Outpus folder path

	extime = clock();

	RunModel();

	std::cout << "Simulation completed" << endl;

	//int pippo;
	//cin >> pippo;

	return 0;
}
#endif

//---------------------------------------------------------------------------
const string Int2Str(const int x)
{
	ostringstream o;
	if (!(o << x)) return "ERROR";
	return o.str();
}
//---------------------------------------------------------------------------
const string Float2Str(const double x)
{
	ostringstream o;
	if (!(o << x)) return "ERROR";
	return o.str();
}
//---------------------------------------------------------------------------
void RunModel(void) {
	std::cout << "Simulation nr. " << para.SimNr << endl;

	// standard deviation for mutational effects for each trait 
	para.mu_std[0] = para.genot_std2[0] / sqrt(2.0 * (double)para.L);
	para.mu_std[1] = para.genot_std2[1] / sqrt(2.0 * (double)para.L);
	para.mu_std[2] = para.genot_std2[2] / sqrt(2.0 * (double)para.L);
	para.mu_std[3] = para.genot_std2[3] / sqrt(2.0 * (double)para.L);

	string name = dirOut + "Sim" + Int2Str(para.SimNr) + "_Para.txt";
	para.outPara(name);

	if (para.out_int > 0) {
		outPop_header();
		if (para.polyEvol || para.dispEvol) outTrait_header();
	}
	if (para.ind_interval > 0) {
		outInds_header();
		outOffs_header();
	}

#if MUTATIONLOAD
	if (para.indRand_interval > 0) outRanInds_header();

	if (para.PopMut_interval > 0) outPopMut_header();
#endif

	//Loop through REPLICATES --------------------------------------------
	for (r = 0; r < para.rep; r++) {
		std::cout << "rep = " << r << "==================" << endl;

#if MUTATIONLOAD
		if (para.out_mutations) outMut_header();
#endif

		//create a new landscape
		landscape(para.x_max, para.y_max, para.prop_suitable);

		//create the environmental gradient
		//envGradient(para.x_max, para.y_max);

		//Initialisation
		initialisation(para.x_max, para.y_max, para.min_seedX, para.min_seedY, para.max_seedX, para.max_seedY);
		//cout << "initialisation ok" << endl;

		mat_trials = false;

		//Loop through GENERATIONS ------------------------------------
		for (g = 0; g < para.gen; g++) {

			if (g % 50 == 0) {
				std::cout << "gen = " << g << endl;
				extime = clock() - extime;
				std::cout << "time = " << (float)extime / CLOCKS_PER_SEC << " sec" << endl;
				extime = clock();
			}

#if STOPEVO
			if (g == para.exp_start - 1) {
				Npop = 0;
				mean_disp = 0.0;

				for (int x = 0; x < para.x_max; x++) {
					for (int y = 0; y < para.y_max; y++) {
						if (pop[x][y] != NULL) {
							if (pop[x][y]->N > 0) {
								mean_disp += pop[x][y]->mFep;
								Npop++;
							}
						}
					}
				}
				mean_disp /= (double)Npop;

				para.dispEvol = 0;
				para.genot_mean[1] = mean_disp;
				cout << para.genot_mean[1] << endl;
			}
#endif
			

			//Gradient shifting
			//if (para.grad_shifts && g >= para.shift_start && g <= para.shift_end) envGradientShift(para.x_max, para.y_max);			

			//Mating
			mating();
			//cout << "mat ok" << endl;

			if (para.postMating) {
				postmating_dispersal(); //Adult post-mating dispersal
				//cout << "disp ok " << endl;
				birth(); //offspring produciton
				//cout << "brth ok " << endl;
			}
			else {
				birth();
				//cout << "brth ok " << endl;
				premating_dispersal(); //natal dispersal
				//cout << "disp ok " << endl;
			}

			//Density-dependent survival
			survival();
			//cout << "surv ok " << endl;

#if MUTATIONLOAD
			if (para.x_max > 1 || para.y_max > 1) {

				if (g > 0 && (g % para.indRand_interval == 0
					|| (g < 502 && g % 10 == 0)
					|| (g > para.exp_start - 1 && g < para.exp_start + 500 && g % 10 == 0))) {
					betweenPop_mating();
					mat_trials = true;
				}
			}
			else {
				if (g > 0 && (g % para.indRand_interval == 0
					|| (g < 502 && g % 10 == 0)
					|| (g > para.exp_start - 1 && g < para.exp_start + 500 && g % 10 == 0))) {
					mat_trials = true;
				}
			}
#endif

			if (Ntot < 2) break;
		}

		//Delete populations and the landscape
		for (int x = 0; x < para.x_max; x++) {
			for (int y = 0; y < para.y_max; y++) {
				if (pop[x][y] != NULL) {
					pop[x][y]->deleteAdults();
					delete pop[x][y]; pop[x][y] = NULL;
				}
				delete land[x][y]; land[x][y] = NULL;
			}
			delete[] pop[x]; pop[x] = NULL;
			delete[] land[x]; land[x] = NULL;
		}
		delete[] pop; pop = NULL;
		delete[] land; land = NULL;
		
		if (muts.is_open()) muts.close();
	}

	if (pops.is_open()) pops.close();
	if (trait.is_open()) trait.close();
	if (inds.is_open()) inds.close();
	if (offs.is_open()) offs.close();
	if (ranInd.is_open()) ranInd.close();
	if (popmut.is_open()) popmut.close();
}
//---------------------------------------------------------------------------
void landscape(int xmax, int ymax, double prop) {

	int x, y, i, ncells;
	std::uniform_int_distribution<> ranx(0, xmax - 1);
	std::uniform_int_distribution<> rany(0, ymax - 1);
	std::uniform_real_distribution<> kdist(para.Kmin, para.Kmax);

	land = new Landscape **[xmax];
	for (int j = 0; j < xmax; j++) {
		land[j] = new Landscape *[ymax];
		for (int jj = 0; jj < ymax; jj++) {
			land[j][jj] = new Landscape();
		}
	}

	pop = new Population **[xmax];
	for (int j = 0; j < xmax; j++) {
		pop[j] = new Population *[ymax];
		//do not intialise the single pop. - do it only when there is a population
		for (int jj = 0; jj < ymax; jj++) pop[j][jj] = NULL;
	}

	if (prop < 1.0) {
		ncells = (int)(xmax * ymax * prop); //nr. of suitable cells 
		//Distribute suitable cells randomly
		i = 0;
		do {
			do {
				x = ranx(rdgen);
				y = rany(rdgen);
			} while (land[x][y]->suitable == 1);
			land[x][y]->suitable = 1;
			if (para.spat_heterog) land[x][y]->local_K = kdist(rdgen); //spatial heterogeneity 
			else land[x][y]->local_K = para.K;
			land[x][y]->local_K0 = new double(land[x][y]->local_K);
			i++;
		} while (i < ncells);
	}
	else {
		for (x = 0; x < xmax; x++) {
			for (y = 0; y < ymax; y++) {
				land[x][y]->suitable = 1;
				if (para.spat_heterog) land[x][y]->local_K = kdist(rdgen); //spatial heterogenetiy
				else land[x][y]->local_K = para.K;
				land[x][y]->local_K0 = new double(land[x][y]->local_K);
			}
		}
	}
}
//---------------------------------------------------------------------------
void initialisation(int xmax, int ymax, int minSx, int minSy, int maxSx, int maxSy) {

	if (para.sexDisp) {
		//Normal distributions for traits initialisation
		std::normal_distribution<> distrPoly(para.genot_mean[0] / (2.0 * (double)para.L), para.genot_std[0] / sqrt(2.0 * (double)para.L));
		std::normal_distribution<> distrEmig(para.genot_mean[1] / (2.0 * (double)para.L), para.genot_std[1] / sqrt(2.0 * (double)para.L));
		std::normal_distribution<> distrEmigMale(para.genot_mean[2] / (2.0 * (double)para.L), para.genot_std[2] / sqrt(2.0 * (double)para.L));
		std::normal_distribution<> distrDist(para.genot_mean[3] / (2.0 * (double)para.L), para.genot_std[2] / sqrt(2.0 * (double)para.L));
		std::normal_distribution<> distrDistMale(para.genot_mean[4] / (2.0 * (double)para.L), para.genot_std[2] / sqrt(2.0 * (double)para.L));

		//Inintialise populations in all the suitable cells in the defined area
		for (int x = minSx; x < maxSx; x++) {
			for (int y = minSy; y < maxSy; y++) {
				if (land[x][y]->suitable) {
					pop[x][y] = new Population(x, y);
					pop[x][y]->initialise_pop(land[x][y]->local_K, k, para, distrEmig, distrDist, distrPoly, s_mild, neutral, position,
						distrEmigMale, distrDistMale);
				}
			}
		}
	}
	else {
		//Normal distributions for traits initialisation
		std::normal_distribution<> distrPoly(para.genot_mean[0] / (2.0 * (double)para.L), para.genot_std[0] / sqrt(2.0 * (double)para.L));
		std::normal_distribution<> distrEmig(para.genot_mean[1] / (2.0 * (double)para.L), para.genot_std[1] / sqrt(2.0 * (double)para.L));
		std::normal_distribution<> distrDist(para.genot_mean[3] / (2.0 * (double)para.L), para.genot_std[2] / sqrt(2.0 * (double)para.L));

		//Inintialise populations in all the suitable cells in the defined area
		for (int x = minSx; x < maxSx; x++) {
			for (int y = minSy; y < maxSy; y++) {
				if (land[x][y]->suitable) {
					pop[x][y] = new Population(x, y);
					pop[x][y]->initialise_pop(land[x][y]->local_K, k, para, distrEmig, distrDist, distrPoly, s_mild, neutral, position);
				}
			}
		}
	}

	cy_max = para.max_seedY;
	cy_min = 0;
}
//---------------------------------------------------------------------------
//// to adapat for calculation of ID and heterosis
//void generateIDdata(void) {
//	int esc, rdnx1, rdny1, rdnx2, rdny2;
//	int mom, dad;
//	Individuals *ind;
//
//	std::uniform_int_distribution<> ran_x(0, para.x_max-1);
//	std::uniform_int_distribution<> ran_y(0, para.y_max-1);
//
//	//generate 10000 offs from random parents across populations
//	for (int i = 0; i < 10000; i++) {
//		
//		//sample mother
//		esc = 0;
//		rdnx1 = -9; rdny1 = -9;
//		mom = -9;
//		do {
//				rdnx1 = ran_x(rdgen);
//				rdny1 = ran_y(rdgen);
//				esc++;
//				
//		} while ((pop[rdnx1][rdny1] == NULL || pop[rdnx1][rdny1]->Nf < 1) && esc < 100);
//
//		if (esc == 100) break;
//
//		if (rdnx1 > -9 && rdny1 > -9) {
//			std::uniform_int_distribution<> rdn_mom(0, pop[rdnx1][rdny1]->Nf - 1);
//			mom = rdn_mom(rdgen);
//		}
//		//sample father
//		esc = 0;
//		rdnx2 = -9; rdny2 = -9;
//		dad = -9;
//		do{
//				rdnx2 = ran_x(rdgen);
//				rdny2 = ran_y(rdgen);
//				esc++;
//				
//		} while ((pop[rdnx2][rdny2] == NULL || pop[rdnx2][rdny2]->Nm < 1) && esc < 100);
//		if (esc == 100) break;
//		
//		if (rdnx2 > -9 && rdny2 > -9) {
//			std::uniform_int_distribution<> rdn_dad(0, pop[rdnx2][rdny2]->Nm - 1);
//			dad = rdn_dad(rdgen);
//		}
//		if (mom > -9 && dad > -9) {
//			//produce offspring
//			ind = new Individuals(para, 0, 0, 0);
//			inheritance(ind, pop[rdnx1][rdny1]->females[mom], pop[rdnx2][rdny2]->males[dad]);
//			ind->outRanInd(1, r, g, &ranInd);
//			delete ind;
//		}
//	}
//}
//---------------------------------------------------------------------------
void env_stochasticity(int xx, int yy)
{
	//temporal environmental stochasticity
	land[xx][yy]->eps = para.ac * land[xx][yy]->eps + t_stoch(rdgen) * std::sqrt(1.0 - para.ac * para.ac);
	land[xx][yy]->local_K = *land[xx][yy]->local_K0 * (1.0 + land[xx][yy]->eps);
	if (land[xx][yy]->local_K < 1.0) {
		land[xx][yy]->local_K = 0.0;
		land[xx][yy]->suitable = false;
	}
	else {
		if (land[xx][yy]->local_K > para.Kmax) land[xx][yy]->local_K = para.Kmax;
		land[xx][yy]->suitable = true;
	}

}
//---------------------------------------------------------------------------
void local_extinction(int xx, int yy) {

		if (extinct(rdgen)) {
			pop[xx][yy]->deleteAdults();
			pop[xx][yy]->Nf = 0;
			pop[xx][yy]->Nm = 0;
			pop[xx][yy]->N = 0;
			pop[xx][yy]->Noffs = 0;
			delete pop[xx][yy];
			pop[xx][yy] = NULL;
		}
}

//---------------------------------------------------------------------------
void mating(void) {
	bool remate;
	int m, dad;
	int a_males; //number of available mates
	double aa, p_surv; //survival probability 
	Individuals* ind;
	vector<int> a_mates; //vector of indeces of available males
	vector<Individuals>::iterator iter;

	//Loop through the landscape
	for (int x = 0; x < para.x_max; x++) {
		for (int y = cy_min; y < cy_max; y++) {
			if (pop[x][y] != NULL) {
				
				if (pop[x][y]->Nf > 0 && pop[x][y]->Nm > 0) {

					//output adult males
					if (para.ind_interval > 0 && g % para.ind_interval == 0) {
						for (int i = 0; i < pop[x][y]->Nm; i++) pop[x][y]->males[i].outInd(para, r, g, &inds);
					}
					
					if (para.matSys > 0) { //polygyny or promiscuity or monogamy
						if (pop[x][y]->Nm > 1) {
							//available males
							a_males = pop[x][y]->Nm;
							for (int i = 0; i < pop[x][y]->Nm; i++) {
								if (pop[x][y]->males[i].reproduce) a_mates.push_back(i); //vector of available mates
								else a_males--;
							}
							if (para.matSys == 1 && a_males > 0) {
								std::uniform_int_distribution<> mat(0, a_males - 1);
								m = mat(rdgen); //single male to sire all females offsring
							}
						}
						else {
							if (pop[x][y]->males[0].reproduce) {
								m = 0;
								a_males = 1;
								a_mates.push_back(0);
							}
							else a_males = 0;
						}
					}
					
					//Loop through females for mating and survival
					for (iter = pop[x][y]->females.begin(); iter != pop[x][y]->females.end(); iter++) {
						p_surv = 1.0; //survival probability

						if (iter->reproduce) {	

							if (para.matSys == 0) {
								if (pop[x][y]->Nm > 1) {
									//available males
									a_males = pop[x][y]->Nm;
									for (int i = 0; i < pop[x][y]->Nm; i++) {
										if (pop[x][y]->males[i].reproduce) a_mates.push_back(i);
										else a_males--;
									}

									if (a_males > 0) {
										//sample the first male
										std::uniform_int_distribution<> mat(0, a_males - 1);
										m = mat(rdgen);
										iter->mates.push_back(a_mates[m]);
										iter->n_mates++;
										a_males--;
										a_mates.erase(a_mates.begin() + m);

#if POLY == 1 											
										//Remating
										if (para.polyEvol) aa = *iter->traits.p_a_mat;
										else aa = para.genot_mean[0];
										std::bernoulli_distribution rem(exp(-aa * (double)iter->n_mates));
										remate = rem(rdgen);
										while (remate && a_males > 0) {
											//sample additional male
											std::uniform_int_distribution<> mat(0, a_males - 1);
											m = mat(rdgen);
											iter->mates.push_back(a_mates[m]);
											iter->n_mates++;
											a_males--;
											a_mates.erase(a_mates.begin() + m);
											//remate?
											std::bernoulli_distribution rem(exp(-aa * (double)iter->n_mates));
											remate = rem(rdgen);
										}
										//determine survival probability based on the cost of realised multiple mating
										if (para.c_f && iter->n_mates > 1)
											p_surv = exp(-pow(1.0 - (double)iter->n_mates, 2.0) / (2.0 * para.omega_f));
										if (unireal(rdgen) > p_surv) iter->alive = false;
#endif
									}
								}
								else {
									if (pop[x][y]->males[0].reproduce) {
										//mate with the only male available
										iter->mates.push_back(0);
										iter->n_mates++;
									}
								}
							}
							else {
								//polygyny
								if (para.matSys == 1 && a_males > 0) {
									iter->n_mates = 1;
									iter->mates.push_back(a_mates[m]);
								}
								//promiscuity
								if (para.matSys == 2) iter->n_mates = a_males;
								//monogamy
								if (para.matSys == 3 && a_males > 0){									
									std::uniform_int_distribution<> mat(0, a_males - 1);
									m = mat(rdgen);
									iter->mates.push_back(a_mates[m]);									
									iter->n_mates = 1;
									a_males--;
									a_mates.erase(a_mates.begin() + m);									
								}
							}
							if (para.matSys == 0 && !a_mates.empty()) a_mates.clear();

#if MUTATIONLOAD
							//produce 10 offs from random mating within-population for ID and heterosis calculations
							if (mat_trials) {
								if ((y > (int)((cy_max - cy_min) / 2) - 4 && y < (int)((cy_max - cy_min) / 2) + 4)) {
									//father sampling distribution
									std::uniform_int_distribution<> sdad(0, pop[x][y]->Nm - 1);

									for (int i = 0; i < 10; i++) {
										ind = new Individuals(para, Bern(rdgen), x, y);
										dad = sdad(rdgen);

										//inheritance
										inheritance(ind, *iter, pop[x][y]->males[dad]);

										//output individual
										ind->outRanInd(0, r, g, &ranInd);
										ind->deleteInd();

										delete ind;
									}
								}
							}
#endif
						} //end if female reproduce
					}//end of female loop
					if (para.matSys > 0 && !a_mates.empty()) a_mates.clear();
				}
			}			
		}
	}

	mat_trials = false;
}
//-------------------------------------------------------------------------------------------------------------
//Post-mating (pre-birth) dispersal
void postmating_dispersal(void) {
	int max_y_disp;
	int maxy;
	int new_x, new_y;
	double x_rand, y_rand;
	double R1, dist, rndAngle;
	double emp; //emigration probabilty 
	double mdist;
	vector<Individuals>::iterator iter;

	std::uniform_real_distribution<> unireal_disp(0.0, 0.999);
	std::uniform_real_distribution<> unireal_dispB(0.0000001, 1.0);

	if (g < para.exp_start) max_y_disp = para.max_seedY;
	else max_y_disp = para.y_max;

	maxy = cy_max;

	for (int x = 0; x < para.x_max; x++) {
		for (int y = cy_min; y < maxy; y++) {
			if (land[x][y]->suitable && pop[x][y] != NULL) {
				if (pop[x][y]->N > 0) {
					//females
					for (iter = pop[x][y]->females.begin(); iter != pop[x][y]->females.end(); iter++) {
						if (iter->alive) {

							if (para.dispersal) {
								if (para.dispEvol == 1 || para.dispEvol == 3)
									emp = *iter->traits.p_ep_f;
								else emp = para.genot_mean[1];
							}
							else emp = 0.0;

							//Disperses?
							if (unireal(rdgen) < emp) {
								if (unireal(rdgen) > para.dispCost) { //survives?

									if (para.nearest) { //nearest-neighbour
										std::uniform_int_distribution<> sample_xy(-1, 1);
										do {
											do {
												new_x = x + sample_xy(rdgen);
												new_y = y + sample_xy(rdgen);
											} while (new_x == x && new_y == y);
										} while (new_x < 0 || new_x >(para.x_max - 1) || new_y < cy_min || new_y >(max_y_disp - 1));

									}
									else {
										//dispersal distance
										if (para.dispEvol == 2 || para.dispEvol == 3) mdist = *iter->traits.p_dist_f;
										else mdist = para.genot_mean[3];
										//sample new location
										x_rand = unireal_disp(rdgen);
										y_rand = unireal_disp(rdgen);
										do {
											do {
												R1 = unireal_dispB(rdgen);
												dist = (-1.0 * mdist) * std::log(R1);
												rndAngle = unireal(rdgen) * 2.0 * PI;
												new_x = (int)(dist * cos(rndAngle) / (double)para.resol + x_rand + x);
												new_y = (int)(dist * sin(rndAngle) / (double)para.resol + y_rand + y);
											} while (new_x == x && new_y == y);
										} while (new_x < 0 || new_x >(para.x_max - 1) || new_y < cy_min || new_y >(max_y_disp - 1));
									}

									//settle or die
									pop[x][y]->Nf--;
									pop[x][y]->N--;
									if (land[new_x][new_y]->suitable) {
										iter->dispersed = true;
										if (pop[new_x][new_y] == NULL) pop[new_x][new_y] = new Population(new_x, new_y);
										pop[new_x][new_y]->tmp_females.push_back(*iter);
										pop[new_x][new_y]->Nf++;
										pop[new_x][new_y]->Nm++;

										if (new_y > cy_max) cy_max = new_y + 1; //update maximum y
									}
									else iter->alive = false;
								}
								else {
									pop[x][y]->Nf--;
									pop[x][y]->Nm--;
									iter->alive = false;
								}
							}
							else { //resident
								pop[x][y]->tmp_females.push_back(*iter);
							}
						}
					}
					//males
					for (iter = pop[x][y]->males.begin(); iter != pop[x][y]->males.end(); iter++) {
						if (iter->alive) {
							if (para.dispersal) {
								if (para.dispEvol == 1 || para.dispEvol == 3) {
									if (para.sexDisp) emp = *iter->traits.p_ep_m;
									else emp = *iter->traits.p_ep_f;
								}
								else {
									if (para.sexDisp) emp = para.genot_mean[2];
									emp = para.genot_mean[1];
								}
							}
							else emp = 0.0;

							//Disperses?
							if (unireal(rdgen) < emp) {
								if (unireal(rdgen) > para.dispCost) { //survives

									if (para.nearest) { //nearest-neighbour
										std::uniform_int_distribution<> sample_xy(-1, 1);
										do {
											do {
												new_x = x + sample_xy(rdgen);
												new_y = y + sample_xy(rdgen);
											} while (new_x == x && new_y == y);
										} while (new_x < 0 || new_x >(para.x_max - 1) || new_y < cy_min || new_y >(max_y_disp - 1));

									}
									else {
										//dispersal distance
										if (para.dispEvol == 2 || para.dispEvol == 3) {
											if (para.sexDisp) mdist = *iter->traits.p_dist_m;
											else mdist = *iter->traits.p_dist_f;
										}
										else {
											if (para.sexDisp) mdist = para.genot_mean[4];
											mdist = para.genot_mean[3];
										}

										//sample new location
										x_rand = unireal_disp(rdgen);
										y_rand = unireal_disp(rdgen);
										do {
											do {
												R1 = unireal_dispB(rdgen);
												dist = (-1.0 * mdist) * std::log(R1);
												rndAngle = unireal(rdgen) * 2.0 * PI;
												new_x = (int)(dist * cos(rndAngle) / (double)para.resol + x_rand + x);
												new_y = (int)(dist * sin(rndAngle) / (double)para.resol + y_rand + y);
											} while (new_x == x && new_y == y);
										} while (new_x < 0 || new_x >(para.x_max - 1) || new_y < cy_min || new_y >(max_y_disp - 1));
									}
									//settle or die
									pop[x][y]->Nm--;
									pop[x][y]->N--;
									if (land[new_x][new_y]->suitable) {
										iter->dispersed = true;
										if (pop[new_x][new_y] == NULL) pop[new_x][new_y] = new Population(new_x, new_y);
										pop[new_x][new_y]->tmp_males.push_back(*iter);
										pop[new_x][new_y]->Nm++;
										pop[new_x][new_y]->N++;

										if (new_y > cy_max) cy_max = new_y + 1; //update maximum y
									}
									else iter->alive = false;
								}
								else {
									pop[x][y]->Nm--;
									pop[x][y]->N--;
									iter->alive = false;
								}
							}
							else { //resident
								pop[x][y]->tmp_males.push_back(*iter);
							}
						}
					}
				}
			}
		}
	}
}

//-------------------------------------------------------------------------------------------------------------
void birth(void) {
	int n_offs;
	int n_mut;
	int dad;
	int n_offs_b; //number of offspring before any selection to calculate frequency of deleterious mutations
	double p_surv;
	Individuals* ind;
	vector<Individuals>::iterator iter;

	//Normal distributions for trait mutations
	std::normal_distribution<> MutPoly(0.0, para.mu_std[0]);
	std::normal_distribution<> MutEm(0.0, para.mu_std[1]);
	std::normal_distribution<> MutDist(0.0, para.mu_std[2]);
	std::normal_distribution<> MutPostEm(0.0, para.mu_std[3]);

	std::poisson_distribution<> n_mildmut(para.Ud);
	std::poisson_distribution<> n_lethal(para.Ul);
	std::poisson_distribution<> n_backmut(para.Ub);


	//Loop through the landscape
	for (int x = 0; x < para.x_max; x++) {
		for (int y = cy_min; y < cy_max; y++) {
		
			if (pop[x][y] != NULL) {

				if (para.postMating == false) { //if there hasn't been post-mating (pre-birth) dispersal
					pop[x][y]->tmp_females = pop[x][y]->females;
					pop[x][y]->tmp_males = pop[x][y]->males;
				}

				/*cout << x<< "\t" << y << endl;
				cout << "shift ok" << endl;*/
				pop[x][y]->Noffs = 0;
				pop[x][y]->Foffs = 0;
				pop[x][y]->Moffs = 0;

				//clear Map of population's deleterious mutations
				if (!pop[x][y]->popMuts.empty()) pop[x][y]->popMuts.clear();
				n_offs_b = 0;

				if (pop[x][y]->Nf > 0 && pop[x][y]->Nm > 0) {

					//Loop through females for mating and survival
					for (iter = pop[x][y]->tmp_females.begin(); iter != pop[x][y]->tmp_females.end(); iter++) {						

						if (iter->alive && iter->reproduce && iter->n_mates > 0) {							
							//number of offspring			
							n_offs = pois(rdgen);
							pop[x][y]->Noffs += n_offs;

							//father sampling distribution (for mating system 0)
							std::uniform_int_distribution<> sdad(0, iter->n_mates - 1);

							if (para.matSys == 1 || para.matSys == 3) dad = iter->mates[0]; //polygyny or monogamy

							//cout << "n mates = " << iter->n_mates << endl;
							//loop through offspring
							for (int i = 0; i < n_offs; i++) {

								ind = new Individuals(para, Bern(rdgen), x, y);

								if (para.matSys == 0 || para.matSys == 2) { //weak polygyny or promiscuity
									if (iter->n_mates > 1) dad = iter->mates[sdad(rdgen)];
									else dad = iter->mates[0];
								}

								/*cout << "dad ok " << dad << endl;
								cout << pop[iter->xnatal][iter->ynatal]->Nm << endl;
								cout << pop[x][y]->Nm << endl;
								cout << land[iter->xnatal][iter->ynatal]->local_K << endl;
								cout << pop[iter->xnatal][iter->ynatal]->males[dad].alive << endl;*/

								//inheritance
								if (para.postMating) inheritance(ind, *iter, pop[iter->xnatal][iter->ynatal]->males[dad]);
								else inheritance(ind, *iter, pop[x][y]->males[dad]);
								
								
								//cout << "inher ok " << endl;

								//neutral mutation
								n_mut = n_neutmut(rdgen);
								if (n_mut > 0) ind->neutral_mutation(n_mut, neutr_position, neutral);
#if MUTATIONLOAD

								////beneficial mutations
								//if (para.ben_prop > 0.0) {
								//	std::poisson_distribution<> n_benmut(para.Ud * para.ben_prop);
								//	n_mut = n_benmut(rdgen);
								//	if (n_mut > 0) ind->benef_mutation(n_mut, para.sb, position);
								//}
								////back mutations
								//n_mut = n_backmut(rdgen);
								//if (n_mut > 0 && ind->chromo.nMut > n_mut) ind->back_mutation(n_mut);

								//mildly deleterius mutations
								n_mut = n_mildmut(rdgen);
								if (n_mut > 0) ind->delet_mutation(n_mut, k, position, s_mild);

								//lethal mutations
								n_mut = n_lethal(rdgen);
								if (n_mut > 0) ind->delet_mutation(n_mut, position, para.sl, para.hl);
#endif
								//mutations to adaptive alleles
								if (g > para.start_dispEvol - 1) adaptiveMutations(ind, MutPoly, MutEm, MutDist, MutPostEm);

								n_offs_b++;

								//cout << "mut ok " << endl;
								if (ind->w > 0.0) {
#if MUTATIONLOAD
									//if mutation load is affecting offspring survival 
									if (para.loadEffect == 0) {
										if (ind->w > 1.0) p_surv = 1.0;
										else p_surv = ind->w;
										std::bernoulli_distribution survOffs(p_surv);
										ind->alive = survOffs(rdgen);
									}

									if (para.loadEffect == -9) {
										p_surv = std::exp(-1.0 * (ind->h / (double)para.nL));
										std::bernoulli_distribution survOffs(p_surv);
										ind->alive = survOffs(rdgen);
									}

									//output offspring
									if (para.ind_interval > 0 && g % para.ind_interval == 0) ind->outOffs(para, r, g, &offs);
#endif
										if (ind->alive) {
											if (ind->sex) {
												pop[x][y]->Jfemales.push_back(*ind);
												pop[x][y]->Foffs++;
											}
											else {
												pop[x][y]->Jmales.push_back(*ind);
												pop[x][y]->Moffs++;
											}
										}
										else {
											ind->deleteInd();
											pop[x][y]->Noffs--;
										}
								}
								else {
									//output offspring
									if (para.ind_interval > 0 && g % para.ind_interval == 0) ind->outOffs(para, r, g, &offs);

									ind->deleteInd();
									pop[x][y]->Noffs--;
								}

								delete ind;

								//cout << "surviv ok " << endl;
											
							}//end new offspring loop												
						} //end if female reproduce
						//output adult female
						if (para.ind_interval > 0 && g % para.ind_interval == 0) iter->outInd(para, r, g, &inds);
					}//end of female loop
#if MUTATIONLOAD
					//output frequency of deleterious mutations
					if (g > para.out_start - 1 && g > 0 && (g % para.PopMut_interval == 0
						|| (g > para.exp_start - 1 && g < para.exp_start + 500 && g % 20 == 0))) {
						pop[x][y]->outMutations(n_offs_b, r, g, &popmut);
					}
#endif
				}
			}
		}
	}
	mat_trials = false;
	//cy_max++;

	//cout << "Start delete" << endl;
	//Loop through the landscape
	for (int x = 0; x < para.x_max; x++) {
		for (int y = cy_min; y < cy_max; y++) {
			if (pop[x][y] != NULL) {
				//all adults die
				pop[x][y]->N = 0;
				pop[x][y]->Nf = 0;
				pop[x][y]->Nm = 0;
				pop[x][y]->deleteAdults();
			}
		}
	}
	//cout << "delete ok " << endl;
}

//---------------------------------------------------------------------------
void inheritance(Individuals *pup, Individuals mom, Individuals dad) {
	int rdn, rdn2, pos, pos2;
	int hom;
	int n_crossovers;
	double cross;
	std::map<double, mutation>::iterator iter, iter2;
	std::set <double> recomSites;
	std::set<double>::iterator itercross;


	//Neutral loci and homozygosity
	rdn = Bern(rdgen);
	rdn2 = Bern(rdgen);

	//SLOWER----------------------------------
	//for (int i = 0; i < para.nL; i++) {
	//	if (rdn) rdn -= recomb(rdgen);
	//	else rdn += recomb(rdgen);
	//	if (rdn2) rdn2 -= recomb(rdgen);
	//	else rdn2 += recomb(rdgen);

	//	pup->markers.push_back(mom.markers[i * 2 + rdn]);
	//	pup->markers.push_back(dad.markers[i * 2 + rdn2]);
	//	if (mom.markers[i * 2 + rdn] == dad.markers[i * 2 + rdn2]) pup->h++;
	//}
	//----------------------------------------

	//FASTER ---------------------------------
	pos = 0; pos2 = 0;
	for (int i = 0; i < para.nL; i++) {
		if (pos == i) {
			if (rdn) rdn = 0;
			else rdn = 1;
			pos += 1 + geom(rdgen);
		}
		if (pos2 == i) {
			if (rdn2) rdn2 = 0;
			else rdn2 = 1;
			pos2 += 1 + geom(rdgen);
		}
		pup->markers.push_back(mom.markers[i * 2 + rdn]);
		pup->markers.push_back(dad.markers[i * 2 + rdn2]);
		if (mom.markers[i * 2 + rdn] == dad.markers[i * 2 + rdn2]) pup->h++;
	}
	//--------------------------------------------

#if ONECHROM
#else
	//Deleterious mutations and adaptive traits are stored on different chromosomes

#if MUTATIONLOAD
	//Recombination of deleterious mutations (see Roze & Rousset 2009, JEB - Appendix 2)
	//Inherit homologue 1 from the mother----------------------------------------------
	if (mom.chromo.nMut > 0) {
		//sample no. of crossovers
		n_crossovers = crossn(rdgen);
		//sample crossover positions
		for (int i = 0; i < n_crossovers; i++) {
			cross = position(rdgen);
			recomSites.insert(cross);
		}
		itercross = recomSites.begin(); //iterator through crossover positions

		//sample starting homologue
		hom = Bern(rdgen);
		iter = mom.chromo.mutations.begin();

		//no mutations before cross-overs positions
		while (n_crossovers > 0 && *itercross < iter->first) {
			itercross++;
			n_crossovers--;
		}

		for (iter = mom.chromo.mutations.begin(); iter != mom.chromo.mutations.end(); iter++) {
			//cross-overs
			while (n_crossovers > 0 && *itercross < iter->first) {
				if (hom == 0) hom++;
				else hom--;
				itercross++;
				n_crossovers--;
			}
			//if mutation is on the right homologue inherit it, otherwise ignore it
			if (iter->second.homol == hom || iter->second.homol == 2) {
				pup->chromo.mutations[iter->first] = iter->second;
				pup->chromo.mutations[iter->first].homol = 0; //inherit first homologue from mom
				pup->chromo.nMut++; 

				//calculate fitness considering the mutation as heterozygote
				pup->w *= (1.0 - iter->second.h * iter->second.s);

				//add mutation to population mutation map
				if (g > para.out_start - 1 && g > 0 && (g % para.PopMut_interval == 0
					|| (g > para.exp_start && g < para.exp_start + 1000 && g % 20 == 0)))
						pop[pup->x][pup->y]->addMutation(iter->first, iter->second.s, iter->second.h);
			}
		}
		if (!recomSites.empty()) recomSites.clear();
	}

	if (dad.chromo.nMut > 0) {
		//sample no. of crossovers
		n_crossovers = crossn(rdgen);
		//sample crossover positions
		for (int i = 0; i < n_crossovers; i++) {
			cross = position(rdgen);
			recomSites.insert(cross);
		}
		itercross = recomSites.begin(); //iterator through crossover positions
										//sample starting homologue
		hom = Bern(rdgen);
		iter = dad.chromo.mutations.begin();

		//no mutations before cross-overs positions
		while (n_crossovers > 0 && *itercross < iter->first) {
			itercross++;
			n_crossovers--;
		}

		for (iter = dad.chromo.mutations.begin(); iter != dad.chromo.mutations.end(); iter++) {
			//crossovers
			while (n_crossovers > 0 && *itercross < iter->first) {
				if (hom == 0) hom++;
				else hom--;
				itercross++;
				n_crossovers--;
			}

			//if mutation is on the right homologue inherit it, otherwise ignore it
			if (iter->second.homol == hom || iter->second.homol == 2) {
				
				iter2 = pup->chromo.mutations.find(iter->first);
				//if mutation is already present --> it is homozygous
				if (iter2 != pup->chromo.mutations.end()) {
					iter2->second.homol = 2; //mutation is homozygote
					pup->chromo.Nho++; 
					//change fitness effect
					pup->w /= (1.0 - iter2->second.h * iter2->second.s);
					pup->w *= (1.0 - iter2->second.s);
				}
				else { //mutation is heterozygote
					pup->chromo.mutations[iter->first] = iter->second;
					pup->chromo.mutations[iter->first].homol = 1;
					pup->chromo.nMut++;
					//fitness effect
					pup->w *= (1.0 - iter->second.h * iter->second.s);
				}
				//add mutation to population mutation map
				if (g > para.out_start - 1 && g > 0 && (g % para.PopMut_interval == 0
					|| (g > para.exp_start && g < para.exp_start + 1000 && g % 20 == 0)))
						pop[pup->x][pup->y]->addMutation(iter->first, iter->second.s, iter->second.h);
			}
		}
		if (!recomSites.empty()) recomSites.clear();
	}

#endif

	//Recombination of the adaptive traits
	if (para.polyEvol || para.dispEvol) {

		if (para.polyEvol) {
			rdn = Bern(rdgen);
			rdn2 = Bern(rdgen);
			for (int i = 0; i < para.L; i++) {
				if (rdn) rdn -= recomb(rdgen);
				else rdn += recomb(rdgen);
				if (rdn2) rdn2 -= recomb(rdgen);
				else rdn2 += recomb(rdgen);
				pup->traits.a_mat[i * 2] = mom.traits.a_mat[i * 2 + rdn];
				pup->traits.a_mat[i * 2 + 1] = dad.traits.a_mat[i * 2 + rdn2];
				*pup->traits.g_a_mat += pup->traits.a_mat[i * 2] + pup->traits.a_mat[i * 2 + 1];
			}
			//phenotype
			*pup->traits.p_a_mat = *pup->traits.g_a_mat;
			if (*pup->traits.p_a_mat < 0.0) *pup->traits.p_a_mat = 0.0;
		}
		if (para.dispEvol > 0) {
			switch (para.dispEvol) {
			case 1:
				rdn = Bern(rdgen);
				rdn2 = Bern(rdgen);
				for (int i = 0; i < para.L; i++) {
					if (rdn) rdn -= recomb(rdgen);
					else rdn += recomb(rdgen);
					if (rdn2) rdn2 -= recomb(rdgen);
					else rdn2 += recomb(rdgen);
					pup->traits.ep_f[i * 2] = mom.traits.ep_f[i * 2 + rdn];
					pup->traits.ep_f[i * 2 + 1] = dad.traits.ep_f[i * 2 + rdn2];
					*pup->traits.g_ep_f += pup->traits.ep_f[i * 2] + pup->traits.ep_f[i * 2 + 1];
				}
				//phenotype
				*pup->traits.p_ep_f = *pup->traits.g_ep_f;
				if (*pup->traits.p_ep_f < 0.0) *pup->traits.p_ep_f = 0.0;
				if (*pup->traits.p_ep_f > 1.0) *pup->traits.p_ep_f = 1.0;
				break;
			case 2:
				rdn = Bern(rdgen);
				rdn2 = Bern(rdgen);
				for (int i = 0; i < para.L; i++) {
					if (rdn) rdn -= recomb(rdgen);
					else rdn += recomb(rdgen);
					if (rdn2) rdn2 -= recomb(rdgen);
					else rdn2 += recomb(rdgen);
					pup->traits.dist_f[i * 2] = mom.traits.dist_f[i * 2 + rdn];
					pup->traits.dist_f[i * 2 + 1] = dad.traits.dist_f[i * 2 + rdn2];
					*pup->traits.g_dist_f = pup->traits.dist_f[i * 2] + pup->traits.dist_f[i * 2 + 1];
				}
				//phenotype
				*pup->traits.p_dist_f = *pup->traits.g_dist_f;
				if (*pup->traits.p_dist_f < 0.0) *pup->traits.p_dist_f = 0.0;
				break;
			case 3:
				rdn = Bern(rdgen);
				rdn2 = Bern(rdgen);
				for (int i = 0; i < para.L; i++) {
					if (rdn) rdn -= recomb(rdgen);
					else rdn += recomb(rdgen);
					if (rdn2) rdn2 -= recomb(rdgen);
					else rdn2 += recomb(rdgen);
					pup->traits.ep_f[i * 2] = mom.traits.ep_f[i * 2 + rdn];
					pup->traits.ep_f[i * 2 + 1] = dad.traits.ep_f[i * 2 + rdn2];
					*pup->traits.g_ep_f += pup->traits.ep_f[i * 2] + pup->traits.ep_f[i * 2 + 1];
				}
				//phenotype
				*pup->traits.p_ep_f = *pup->traits.g_ep_f;
				if (*pup->traits.p_ep_f < 0.0) *pup->traits.p_ep_f = 0.0;
				if (*pup->traits.p_ep_f > 1.0) *pup->traits.p_ep_f = 1.0;

				rdn = Bern(rdgen);
				rdn2 = Bern(rdgen);
				for (int i = 0; i < para.L; i++) {
					if (rdn) rdn -= recomb(rdgen);
					else rdn += recomb(rdgen);
					if (rdn2) rdn2 -= recomb(rdgen);
					else rdn2 += recomb(rdgen);
					pup->traits.dist_f[i * 2] = mom.traits.dist_f[i * 2 + rdn];
					pup->traits.dist_f[i * 2 + 1] = dad.traits.dist_f[i * 2 + rdn2];
					*pup->traits.g_dist_f = pup->traits.dist_f[i * 2] + pup->traits.dist_f[i * 2 + 1];
				}
				//phenotype
				*pup->traits.p_dist_f = *pup->traits.g_dist_f;
				if (*pup->traits.p_dist_f < 0.0) *pup->traits.p_dist_f = 0.0;

				break;
			}
			if (para.sexDisp) {
				switch (para.dispEvol) {
				case 1:
					rdn = Bern(rdgen);
					rdn2 = Bern(rdgen);
					for (int i = 0; i < para.L; i++) {
						if (rdn) rdn -= recomb(rdgen);
						else rdn += recomb(rdgen);
						if (rdn2) rdn2 -= recomb(rdgen);
						else rdn2 += recomb(rdgen);
						pup->traits.ep_m[i * 2] = mom.traits.ep_m[i * 2 + rdn];
						pup->traits.ep_m[i * 2 + 1] = dad.traits.ep_m[i * 2 + rdn2];
						*pup->traits.g_ep_m += pup->traits.ep_m[i * 2] + pup->traits.ep_m[i * 2 + 1];
					}
					//phenotype
					*pup->traits.p_ep_m = *pup->traits.g_ep_m;
					if (*pup->traits.p_ep_m < 0.0) *pup->traits.p_ep_m = 0.0;
					if (*pup->traits.p_ep_m > 1.0) *pup->traits.p_ep_m = 1.0;
					break;
				case 2:
					rdn = Bern(rdgen);
					rdn2 = Bern(rdgen);
					for (int i = 0; i < para.L; i++) {
						if (rdn) rdn -= recomb(rdgen);
						else rdn += recomb(rdgen);
						if (rdn2) rdn2 -= recomb(rdgen);
						else rdn2 += recomb(rdgen);
						pup->traits.dist_m[i * 2] = mom.traits.dist_m[i * 2 + rdn];
						pup->traits.dist_m[i * 2 + 1] = dad.traits.dist_m[i * 2 + rdn2];
						*pup->traits.g_dist_m = pup->traits.dist_m[i * 2] + pup->traits.dist_m[i * 2 + 1];
					}
					//phenotype
					*pup->traits.p_dist_m = *pup->traits.g_dist_m;
					if (*pup->traits.p_dist_m < 0.0) *pup->traits.p_dist_m = 0.0;
					break;
				case 3:
					rdn = Bern(rdgen);
					rdn2 = Bern(rdgen);
					for (int i = 0; i < para.L; i++) {
						if (rdn) rdn -= recomb(rdgen);
						else rdn += recomb(rdgen);
						if (rdn2) rdn2 -= recomb(rdgen);
						else rdn2 += recomb(rdgen);
						pup->traits.ep_m[i * 2] = mom.traits.ep_m[i * 2 + rdn];
						pup->traits.ep_m[i * 2 + 1] = dad.traits.ep_m[i * 2 + rdn2];
						*pup->traits.g_ep_m += pup->traits.ep_m[i * 2] + pup->traits.ep_m[i * 2 + 1];
					}
					//phenotype
					*pup->traits.p_ep_m = *pup->traits.g_ep_m;
					if (*pup->traits.p_ep_m < 0.0) *pup->traits.p_ep_m = 0.0;
					if (*pup->traits.p_ep_m > 1.0) *pup->traits.p_ep_m = 1.0;

					rdn = Bern(rdgen);
					rdn2 = Bern(rdgen);
					for (int i = 0; i < para.L; i++) {
						if (rdn) rdn -= recomb(rdgen);
						else rdn += recomb(rdgen);
						if (rdn2) rdn2 -= recomb(rdgen);
						else rdn2 += recomb(rdgen);
						pup->traits.dist_m[i * 2] = mom.traits.dist_m[i * 2 + rdn];
						pup->traits.dist_m[i * 2 + 1] = dad.traits.dist_m[i * 2 + rdn2];
						*pup->traits.g_dist_m = pup->traits.dist_m[i * 2] + pup->traits.dist_m[i * 2 + 1];
					}
					//phenotype
					*pup->traits.p_dist_m = *pup->traits.g_dist_m;
					if (*pup->traits.p_dist_m < 0.0) *pup->traits.p_dist_m = 0.0;
					break;
				}
			}
		}
	}

#endif
}
//---------------------------------------------------------------------------
void adaptiveMutations(Individuals *pup, std::normal_distribution<> MutPoly, std::normal_distribution<> MutEm, 
	std::normal_distribution<> MutDist, std::normal_distribution<> MutPostEm) {
	int nmut;
	int allele;

	if (para.polyEvol || para.dispEvol) {
		if (para.polyEvol) {
			nmut = n_amut(rdgen);
			if (nmut > 0) {
				for (int i = 0; i < nmut; i++) {
					//sample allele
					allele = uni_loci(rdgen);
					pup->traits_mutation(0, allele, MutPoly);
				}
			}
		}
		if (para.dispEvol > 0) {
			switch (para.dispEvol) {
			case 1:
				nmut = n_amut(rdgen);
				if (nmut > 0) {
					for (int i = 0; i < nmut; i++) {
						allele = uni_loci(rdgen);
						
						pup->traits_mutation(1, allele, MutEm);
					}
				}
				break;
			case 2:
				nmut = n_amut(rdgen);
				if (nmut > 0) {
					for (int i = 0; i < nmut; i++) {
						allele = uni_loci(rdgen);
						pup->traits_mutation(2, allele, MutDist);
					}
				}
				break;
			case 3:
				nmut = n_amut(rdgen);
				if (nmut > 0) {
					for (int i = 0; i < nmut; i++) {
						allele = uni_loci(rdgen);
						pup->traits_mutation(1, allele, MutEm);
					}
				}
				nmut = n_amut(rdgen);
				if (nmut > 0) {
					for (int i = 0; i < nmut; i++) {
						allele = uni_loci(rdgen);
						pup->traits_mutation(2, allele, MutDist);
					}
				}
				break;
			}
			if (para.sexDisp) {
				switch (para.dispEvol) {
				case 1:
					nmut = n_amut(rdgen);
					if (nmut > 0) {
						for (int i = 0; i < nmut; i++) {
							allele = uni_loci(rdgen);
							pup->traits_mutation(3, allele, MutEm);
						}
					}
					break;
				case 2:
					nmut = n_amut(rdgen);
					if (nmut > 0) {
						for (int i = 0; i < nmut; i++) {
							allele = uni_loci(rdgen);
							pup->traits_mutation(4, allele, MutDist);
						}
					}
					break;
				case 3:
					nmut = n_amut(rdgen);
					if (nmut > 0) {
						for (int i = 0; i < nmut; i++) {
							allele = uni_loci(rdgen);
							pup->traits_mutation(3, allele, MutEm);
						}
					}
					nmut = n_amut(rdgen);
					if (nmut > 0) {
						for (int i = 0; i < nmut; i++) {
							allele = uni_loci(rdgen);
							pup->traits_mutation(4, allele, MutDist);
						}
					}
					break;
				}
			}
		}
	}
}
//--------------------------------------------------------------------------
//natal dispersal ----------------------------------------------------------
void premating_dispersal(void) {
	int max_y_disp;
	int maxy;
	int new_x, new_y;
	double x_rand, y_rand;
	double R1, dist, rndAngle;
	double emp; //emigration probabilty 
	double mdist;
	vector<Individuals>::iterator iter;

	std::uniform_real_distribution<> unireal_disp(0.0, 0.999);
	std::uniform_real_distribution<> unireal_dispB(0.0000001, 1.0);

	if (g < para.exp_start) max_y_disp = para.max_seedY;
	else max_y_disp = para.y_max;

	maxy = cy_max;

	for (int x = 0; x < para.x_max; x++) {
		for (int y = cy_min; y < maxy; y++) {
			if (land[x][y]->suitable && pop[x][y] != NULL) {

				if (pop[x][y]->Noffs > 0) {
					//females
					for (iter = pop[x][y]->Jfemales.begin(); iter != pop[x][y]->Jfemales.end(); iter++) {
						if (para.dispersal) {
							if (para.dispEvol == 1 || para.dispEvol == 3)
								emp = *iter->traits.p_ep_f;
							else emp = para.genot_mean[1];
						}
						else emp = 0.0;

						//Disperses?
						if (unireal(rdgen) < emp){
							if (unireal(rdgen) > para.dispCost) { //survives?
								
								if (para.nearest) { //nearest-neighbour
									std::uniform_int_distribution<> sample_xy(-1, 1);
									do {
										do {
											new_x = x + sample_xy(rdgen);
											new_y = y + sample_xy(rdgen);
										} while (new_x == x && new_y == y);
									} while (new_x < 0 || new_x >(para.x_max - 1) || new_y < cy_min || new_y >(max_y_disp - 1));

								}
								else {
									//dispersal distance
									if (para.dispEvol == 2 || para.dispEvol == 3) mdist = *iter->traits.p_dist_f;
									else mdist = para.genot_mean[3];
									//sample new location
									x_rand = unireal_disp(rdgen);
									y_rand = unireal_disp(rdgen);
									do {
										do {
											R1 = unireal_dispB(rdgen);
											dist = (-1.0 * mdist) * std::log(R1);
											rndAngle = unireal(rdgen) * 2.0 * PI;
											new_x = (int)(dist * cos(rndAngle) / (double)para.resol + x_rand + x);
											new_y = (int)(dist * sin(rndAngle) / (double)para.resol + y_rand + y);
										} while (new_x == x && new_y == y);
									} while (new_x < 0 || new_x >(para.x_max - 1) || new_y < cy_min || new_y >(max_y_disp - 1));
								}

								//settle or die
								pop[x][y]->Foffs--;
								pop[x][y]->Noffs--;
								if (land[new_x][new_y]->suitable) {
									iter->dispersed = true;
									if (pop[new_x][new_y] == NULL) pop[new_x][new_y] = new Population(new_x, new_y);
									pop[new_x][new_y]->tmp_females.push_back(*iter);
									pop[new_x][new_y]->Foffs++;
									pop[new_x][new_y]->Noffs++;

									if (new_y > cy_max) cy_max = new_y + 1; //update maximum y
								}
								else iter->deleteInd();
							}
							else {
								pop[x][y]->Foffs--;
								pop[x][y]->Noffs--;
								iter->deleteInd();
							}
						}
						else { //resident
							pop[x][y]->tmp_females.push_back(*iter);
						}
					}
					//males
					for (iter = pop[x][y]->Jmales.begin(); iter != pop[x][y]->Jmales.end(); iter++) {
						if (para.dispersal) {
							if (para.dispEvol == 1 || para.dispEvol == 3) {
								if (para.sexDisp) emp = *iter->traits.p_ep_m;
								else emp = *iter->traits.p_ep_f;
							}
							else {
								if (para.sexDisp) emp = para.genot_mean[2];
								emp = para.genot_mean[1];
							}
						}
						else emp = 0.0;
				
						//Disperses?
						if (unireal(rdgen) < emp) {
							if (unireal(rdgen) > para.dispCost) { //survives

								if (para.nearest) { //nearest-neighbour
									std::uniform_int_distribution<> sample_xy(-1, 1);
									do {
										do {
											new_x = x + sample_xy(rdgen);
											new_y = y + sample_xy(rdgen);
										} while (new_x == x && new_y == y);
									} while (new_x < 0 || new_x >(para.x_max - 1) || new_y < cy_min || new_y >(max_y_disp - 1));

								}
								else {
									//dispersal distance
									if (para.dispEvol == 2 || para.dispEvol == 3) {
										if (para.sexDisp) mdist = *iter->traits.p_dist_m;
										else mdist = *iter->traits.p_dist_f;
									}
									else {
										if (para.sexDisp) mdist = para.genot_mean[4];
										mdist = para.genot_mean[3];
									}

									//sample new location
									x_rand = unireal_disp(rdgen);
									y_rand = unireal_disp(rdgen);
									do {
										do {
											R1 = unireal_dispB(rdgen);
											dist = (-1.0 * mdist) * std::log(R1);
											rndAngle = unireal(rdgen) * 2.0 * PI;
											new_x = (int)(dist * cos(rndAngle) / (double)para.resol + x_rand + x);
											new_y = (int)(dist * sin(rndAngle) / (double)para.resol + y_rand + y);
										} while (new_x == x && new_y == y);
									} while (new_x < 0 || new_x >(para.x_max - 1) || new_y < cy_min || new_y >(max_y_disp - 1));
								}
								//settle or die
								pop[x][y]->Moffs--;
								pop[x][y]->Noffs--;
								if (land[new_x][new_y]->suitable) {
									iter->dispersed = true;
									if (pop[new_x][new_y] == NULL) pop[new_x][new_y] = new Population(new_x, new_y);
									pop[new_x][new_y]->tmp_males.push_back(*iter);
									pop[new_x][new_y]->Moffs++;
									pop[new_x][new_y]->Noffs++;

									if (new_y > cy_max) cy_max = new_y + 1; //update maximum y
								}
								else iter->deleteInd();
							}
							else {
								pop[x][y]->Moffs--;
								pop[x][y]->Noffs--;
								iter->deleteInd();
							}
						}
						else { //resident
							pop[x][y]->tmp_males.push_back(*iter);
						}
					}
				}
				pop[x][y]->Jfemales.clear();
				pop[x][y]->Jmales.clear();
			}
		}
	}
}
//--------------------------------------------------------------------------
//Density-dependent survival
void survival(void) {
	int maxy;
	double ps, pr;
	vector<Individuals>::iterator iter;

	Ntot = 0;

	maxy = cy_max;
	//delete populations that fall off the minimum y
	if (cy_min < (cy_max - para.cutoff)) {
		for (int x = 0; x < para.x_max; x++) {
			for (int y = cy_min; y < (cy_max - para.cutoff); y++) {
				if (land[x][y]->suitable && pop[x][y] != NULL) {
					if (para.postMating) {
						for (iter = pop[x][y]->Jfemales.begin(); iter != pop[x][y]->Jfemales.end(); iter++) iter->deleteInd();
						for (iter = pop[x][y]->Jmales.begin(); iter != pop[x][y]->Jmales.end(); iter++) iter->deleteInd();
						pop[x][y]->Jfemales.clear();
						pop[x][y]->Jmales.clear();
					}
					else {
						for (iter = pop[x][y]->tmp_females.begin(); iter != pop[x][y]->tmp_females.end(); iter++) iter->deleteInd();
						for (iter = pop[x][y]->tmp_males.begin(); iter != pop[x][y]->tmp_males.end(); iter++) iter->deleteInd();
						pop[x][y]->tmp_females.clear();
						pop[x][y]->tmp_males.clear();
					}
					delete pop[x][y];
					pop[x][y] = NULL;
				}
			}
		}
		cy_min = cy_max - para.cutoff;
	}

	for (int x = 0; x < para.x_max; x++) {
		for (int y = cy_min; y < maxy; y++) {
			
			if (land[x][y]->suitable && pop[x][y] != NULL) {
				if (pop[x][y]->Noffs > 0) {
						
					if (para.postMating) {
						pop[x][y]->tmp_females = pop[x][y]->Jfemales;
						pop[x][y]->tmp_males = pop[x][y]->Jmales;
						pop[x][y]->Jfemales.clear();
						pop[x][y]->Jmales.clear();
					}

					pop[x][y]->set2zero(); //set population stats to zero

					ps = std::fmin(land[x][y]->local_K / (double)pop[x][y]->Noffs, 1.0);
					std::bernoulli_distribution survive(ps);

					for (iter = pop[x][y]->tmp_females.begin(); iter != pop[x][y]->tmp_females.end(); iter++) {
						if (survive(rdgen)) {
							//Mutation load affect probability of reproducing
							if (para.loadEffect == 1) {
								if (iter->w > 1.0) pr = 1.0;
								else pr = iter->w;
								std::bernoulli_distribution repr(pr);
								iter->reproduce = repr(rdgen);
							}
							pop[x][y]->females.push_back(*iter);
							pop[x][y]->Nf++;
							pop[x][y]->N++;
							pop[x][y]->computeSums(para, *iter);
#if MUTATIONLOAD
								
							//output mutations
							if (para.out_mutations && (g % 5000 == 0 
								|| (g > para.exp_start - 500 && g < para.exp_start + 2000 && g % 100 == 0)))
								iter->outMut(g, &muts);
#endif
						}
						else iter->deleteInd();
					}
					for (iter = pop[x][y]->tmp_males.begin(); iter != pop[x][y]->tmp_males.end(); iter++) {
						if (survive(rdgen)) {
							//Mutation load affect probability of reproducing
							if (para.loadEffect == 1) {
								if (iter->w > 1.0) pr = 1.0;
								else pr = iter->w;
								std::bernoulli_distribution repr(pr);
								iter->reproduce = repr(rdgen);
							}
							pop[x][y]->males.push_back(*iter);
							pop[x][y]->Nm++;
							pop[x][y]->N++;
							pop[x][y]->computeSums(para, *iter);

#if MUTATIONLOAD
							//output mutations
							if (para.out_mutations && (g % 5000 == 0
								|| (g > para.exp_start - 500 && g < para.exp_start + 2000 && g % 100 == 0)))
								iter->outMut(g, &muts);
#endif
						}
						else iter->deleteInd();
					}

					if (pop[x][y]->N > 0) {

						if (para.p_cell_ext > 0.0 && extinct(rdgen)) { //Local extinctions
							//all adults die
							pop[x][y]->N = 0;
							pop[x][y]->Nf = 0;
							pop[x][y]->Nm = 0;
							pop[x][y]->deleteAdults();
							delete pop[x][y];
							pop[x][y] = NULL;
						}
						else { //populationn survives
							pop[x][y]->computeStats(para);

							//output populations and traits
							if (g > para.out_start - 1 && (g % para.out_int == 0
								|| (g > para.exp_start - 100 && g < para.exp_start + 1500 && g % 10 == 0))) {
								pop[x][y]->outPop(r, g, &pops);
								if ((para.polyEvol || para.dispEvol) && g > para.start_dispEvol - 1)
									pop[x][y]->outTrait(para, r, g, &trait);
							}
							pop[x][y]->age++;

							pop[x][y]->Noffs = 0;
							pop[x][y]->Foffs = 0;
							pop[x][y]->Moffs = 0;
							pop[x][y]->tmp_females.clear();
							pop[x][y]->tmp_males.clear();

							Ntot += pop[x][y]->N;
						}
					}
					else {
						delete pop[x][y];
						pop[x][y] = NULL;
					}
				}
				else {
					delete pop[x][y];
					pop[x][y] = NULL;
				}			
			}
			//temporal environmental stochasticity
			if (para.env_stoch) env_stochasticity(x, y);
			if (land[x][y]->local_K == 0.0 && pop[x][y] != NULL) {
				//all adults die
				pop[x][y]->N = 0;
				pop[x][y]->Nf = 0;
				pop[x][y]->Nm = 0;
				pop[x][y]->deleteAdults();
				delete pop[x][y];
				pop[x][y] = NULL;
			}
		}
	}
}
//---------------------------------------------------------------------------
void betweenPop_mating(void) {
	int maxy = cy_max;
	int xx, yy;
	int dad;
	Individuals *ind;

	vector<Individuals>::iterator iter;

	std::uniform_int_distribution<> sample_x(0, para.x_max - 1);
	std::uniform_int_distribution<> sample_y(cy_min, cy_max - 1);
	
	//Loop through the landscape
	for (int x = 0; x < para.x_max; x++) {
		for (int y = cy_min; y < maxy; y++) {
			//if population from the core or the front 5 rows
			//if (y > (cy_max - 10) || (y > (int)((cy_max - cy_min) / 2) - 3 && y < (int)((cy_max - cy_min) / 2) + 3)) {
			//---only core
			if ((y > (int)((cy_max - cy_min) / 2) - 4 && y < (int)((cy_max - cy_min) / 2) + 4)) {
				if (pop[x][y] != NULL) {
					if (pop[x][y]->Nf > 0) {
						//Loop through females and generate 10 offspring each with males from other random populations
						for (iter = pop[x][y]->females.begin(); iter != pop[x][y]->females.end(); iter++) {

							for (int i = 0; i < 10; i++) {
								do {
									xx = sample_x(rdgen);
									yy = sample_y(rdgen);
								} while (pop[xx][yy] == NULL || pop[xx][yy]->Nm < 1 || (xx == x && yy == y));
								std::uniform_int_distribution<> sample_mate(0, pop[xx][yy]->Nm - 1);
								dad = sample_mate(rdgen);

								ind = new Individuals(para, Bern(rdgen), x, y);
								//inheritance
								inheritance(ind, *iter, pop[xx][yy]->males[dad]);

								//output individual
								ind->outRanInd(1, r, g+1, &ranInd);

								ind->deleteInd();

								delete ind;
							}
						}
					}
				}
			}
		}
	}

}
//---------------------------------------------------------------------------

void outPop_header(void) {
	string name;

	name = dirOut + "Sim" + Int2Str(para.SimNr) + "_Pops.txt";
	pops.open(name.c_str());
	pops << "rep\tgen\tx\ty\tage\tNf\tNm";
#if MUTATIONLOAD
	pops << "\tW\tWstd\tHomoz\tHstd\tDelMut\tDelMutSd";
#endif
	pops << endl;
}
//----------------------------------------------------------------------------------------
void outTrait_header(void) {
	string name;

	name = dirOut + "Sim" + Int2Str(para.SimNr) + "_Traits.txt";
	trait.open(name.c_str());

	trait << "rep\tgen\tx\ty\tmMates\tsdMates";
	if (para.polyEvol) trait << "\tA\tAstd";
	if (para.dispEvol > 0) {
		switch (para.dispEvol) {
		case 1:
			trait << "\tfep\tfepStd\tgfep\tgfepStd";
			break;
		case 2:
			trait << "\tfdist\tfdistStd\tgfdist\tgfdistStd";
			break;
		case 3:
			trait << "\tfep\tfepStd\tgfep\tgfepStd";
			trait << "\tfdist\tfdistStd\tgfdist\tgfdistStd";
			break;
		}
		if (para.sexDisp) {
			switch (para.dispEvol) {
			case 1:
				trait << "\tmep\tmepStd\tgmep\tgmepStd";
				break;
			case 2:
				trait << "\tmdist\tmdistStd\tgmdist\tgmdistStd";
				break;
			case 3:
				trait << "\tmep\tmepStd\tgmep\tgmepStd";
				trait << "\tmdist\tmdistStd\tgmdist\tgmdistStd";
				break;
			}
		}
	}
	trait << endl;
}
//----------------------------------------------------------------------------------------
void outMut_header(void) {
	string name;

	name = dirOut + "Sim" + Int2Str(para.SimNr) + "_Muts_rep" + Int2Str(r) + ".txt";
	muts.open(name.c_str());
	muts << "gen\ts\th\thomo" << endl;
}
//----------------------------------------------------------------------------------------
void outInds_header(void) {
	string name;

	name = dirOut + "Sim" + Int2Str(para.SimNr) + "_Inds.txt";
	inds.open(name.c_str());
	inds << "rep\tgen\tx\ty\tsex\talive\treproduce\tdispersed\tnmates\tW\thomoz";
#if MUTATIONLOAD
	inds << "\tNmut\tNho";
#endif
	if (para.polyEvol) inds << "\tgA\tpA";
	if (para.dispEvol > 0) {
		switch (para.dispEvol){
		case 1:
			inds << "\tgEp\tpEp";
			break;
		case 2:
			inds << "\tgDist\tpDist";
			break;
		case 3:
			inds << "\tgEp\tpEp\tgDist\tpDist";
			break;
		}
	}
	inds << endl;
}
//----------------------------------------------------------------------------------------
void outOffs_header(void) {
	string name;

	name = dirOut + "Sim" + Int2Str(para.SimNr) + "_Offs.txt";
	offs.open(name.c_str());
	offs << "rep\tgen\tx\ty\tsex\talive\tW\thomoz";
#if MUTATIONLOAD
	offs << "\tNmut\tNho";
#endif
	offs << endl;
}
//----------------------------------------------------------------------------------------
void outRanInds_header(void)
{
	string name;

	name = dirOut + "Sim" + Int2Str(para.SimNr) + "_RanInds.txt";
	ranInd.open(name.c_str());
	//"out" indicates whether it is a "between population" off (1) or a "within population" off (0)
	ranInd << "rep\tgen\tx\ty\tout\tW\thomoz\tNmut\tNho";
	ranInd << endl;
}
//----------------------------------------------------------------------------------------
void outPopMut_header(void)
{
	string name;
	name = dirOut + "Sim" + Int2Str(para.SimNr) + "_PopMut.txt";
	popmut.open(name.c_str());

	popmut << "rep\tgen\tx\ty\ts\th\tfreq" << endl;
}
