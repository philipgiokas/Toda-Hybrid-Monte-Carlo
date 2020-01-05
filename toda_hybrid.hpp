#ifndef TODA_HYBRID_DEF_3
#define TODA_HYBRID_DEF_3


#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <fstream>
#include <cassert>
#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>



using namespace std;


class toda_hybrid_3
{
	public:
                string type;
		vector<vector<vector<long double>>> phi;
		vector<vector<vector<long double>>> phi_old;
		vector<vector<vector<long double>>> pi;
                vector<vector<long double>> bnd_dat;
                vector<vector<long double>> wall_corr_buff;
                vector<long double> wall_corr;
		gsl_rng* rlx = gsl_rng_alloc(gsl_rng_ranlxd2);
                gsl_rng* r48 = gsl_rng_alloc(gsl_rng_rand48);
		int fit_ln ,latt_ln, nmd ,iters, iters_log2, num_buffs_per_bin, 
                num_buffs_per_bin_log2, num_bins, count, attempts, num_smps_per_buff, num_smps_per_buff_log_2, num_buffs;
		long double beta, m, tau;
		long double beta_alpha[3][2], c[3], mn_fld[2] = {0.0, 0.0}; 
                toda_hybrid_3(string, int,long double, long double, long double, int, int, int);
                long double sq(long double);
                long double hamil();
                bool update();
                //void thermalise(int);
                void phi_to_old();
                void old_to_phi();
                void pi_init();
                void pi_hlf_stp();
                void pi_fll_stp();
                void phi_fll_stp();
                void smp_wall_corr();
                void smp_wall_corr_buff(int);
                void calc_wall_corr_frm_buff();
                void crt_btstrp_wall_corr(double* y);
                void calc_wall_corr_frm_bnd_dat(vector<long double> &);
                void calc_err_at_bn_dpth(int, int, double[]);
                void calc_err(double[]);
                void param_set();
                void set_beta();
                void set_tau();
                void smp_mn_fld();
                void prnt_mn_fld();
                void prnt_results();
                void crt_bnd_dat(int);
                ~toda_hybrid_3();
};







#endif