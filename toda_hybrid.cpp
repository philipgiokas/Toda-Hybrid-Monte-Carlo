#include "toda_hybrid.hpp"
#include "mass_curve_fit.hpp"

toda_hybrid_3::toda_hybrid_3(string type, int latt_ln, long double beta, 
	   long double m, long double tau, int nmd, int iters_log2, int fit_ln)
{
	this -> latt_ln = latt_ln;
	this -> m = m;
	this -> beta = beta;
	this -> tau = tau;
	this -> nmd = nmd;
	this -> iters_log2 = iters_log2;
	this -> fit_ln = fit_ln;
	this -> type = type;

	num_smps_per_buff_log_2 = 4;
	num_smps_per_buff = 1 << num_smps_per_buff_log_2;
    iters = 1 << iters_log2;
    num_buffs = iters / num_smps_per_buff;

    param_set();

	phi.resize(2, vector<vector<long double>> (latt_ln,vector<long double>(latt_ln,0.0)));

	phi_old.resize(2, vector<vector<long double>>(latt_ln, vector<long double>(latt_ln, 0.0)));

	pi.resize(2, vector<vector<long double>>(latt_ln,vector<long double>(latt_ln, 0.0)));

	wall_corr_buff.resize(3 * fit_ln, vector<long double>(num_buffs, 0.0));

    wall_corr.resize(3 * fit_ln , 0.0);


	for(int i = 0 ; i < latt_ln ; i++)
	{
		for(int j = 0 ; j < latt_ln ; j++)
		{
			phi[0][i][j] = gsl_ran_gaussian(rlx, 1.0);
			phi[1][i][j] = gsl_ran_gaussian(rlx, 1.0);
		}
	}


	
}

void toda_hybrid_3::param_set()
{

    assert(type == "g" || type == "d" ||
           type == "a" || type == "ai"||
           type == "c" || type == "ci");
    if(type == "g")
    {
    	beta_alpha[0][0] = beta * sqrt(2.0);
		beta_alpha[0][1] = 0.0;
		beta_alpha[1][0] = - beta / sqrt(2.0);
		beta_alpha[1][1] = - beta * sqrt(1.5);
		beta_alpha[2][0] = - beta / sqrt(2.0);
		beta_alpha[2][1] = beta / sqrt(6.0);

		c[0] = 2.0 * sq(m / beta);
		c[1] = sq(m / beta);
		c[2] = 3.0 * sq(m / beta);
    }
    else if (type == "d")
    {
    	beta_alpha[0][0] = beta * sqrt(2.0);
		beta_alpha[0][1] = 0.0;
		beta_alpha[1][0] = - beta / sqrt(2.0);
		beta_alpha[1][1] = - beta * sqrt(1.5);
		beta_alpha[2][0] = - 3.0 * beta / sqrt(2.0);
		beta_alpha[2][1] = beta * sqrt(1.5);
	  
	    c[0] = 2.0 * sq(m / beta);
	    c[1] = sq(m / beta);
	    c[2] = sq(m / beta);
    }
	else if (type == "a")
    {
    	beta_alpha[0][0] = -beta;
		beta_alpha[0][1] =  beta;
		beta_alpha[1][0] =  beta;
		beta_alpha[1][1] =  0.0;
		beta_alpha[2][0] =  -0.5 * beta;
		beta_alpha[2][1] =  -0.5 * beta;
	  
	    c[0] = sq(m / beta);
	    c[1] = 2.0 * sq(m / beta);
	    c[2] = 2.0 * sq(m / beta);
    }
    else if (type == "ai")
    {
    	beta_alpha[0][0] = -beta;
		beta_alpha[0][1] =  beta;
		beta_alpha[1][0] =  2.0 * beta;
		beta_alpha[1][1] =  0.0;
		beta_alpha[2][0] =  -2.0 * beta;
		beta_alpha[2][1] =  -2.0 * beta;
	  
	    c[0] = sq(m / beta);
	    c[1] = sq(m / beta);
	    c[2] = 0.5 * sq(m / beta);
    }
    else if (type == "c")
    {
    	beta_alpha[0][0] = -beta;
		beta_alpha[0][1] =  beta;
		beta_alpha[1][0] =   beta;
		beta_alpha[1][1] =  0.0;
		beta_alpha[2][0] =  - beta;
		beta_alpha[2][1] =  - beta;
	  
	    c[0] = sq(m / beta);
	    c[1] = 2.0 * sq(m / beta);
	    c[2] = sq(m / beta);
    }
    else if (type == "ci")
    {
    	beta_alpha[0][0] = -beta;
		beta_alpha[0][1] =  beta;
		beta_alpha[1][0] =  2.0 * beta;
		beta_alpha[1][1] =  0.0;
		beta_alpha[2][0] =  -beta;
		beta_alpha[2][1] =  -beta;
	  
	    c[0] = sq(m / beta);
	    c[1] = sq(m / beta);
	    c[2] = sq(m / beta);
    }
}

void toda_hybrid_3::pi_init()
{
	for(int i = 0 ; i < latt_ln ; i++)
	{
		for(int j = 0 ; j < latt_ln ; j++)
		{
			pi[0][i][j] = gsl_ran_gaussian(rlx, 1.0);
			pi[1][i][j] = gsl_ran_gaussian(rlx, 1.0);
		}
	}
}

long double toda_hybrid_3::sq(long double in)
{
	return in * in;
}

void toda_hybrid_3::phi_to_old()
{
	for(int i = 0 ; i < latt_ln ; i++)
	{
		for(int j = 0 ; j < latt_ln ; j++)
		{
			phi[0][i][j] = phi_old[0][i][j];
			phi[1][i][j] = phi_old[1][i][j];
		}
	}
}

void toda_hybrid_3::old_to_phi()
{
	for(int i = 0 ; i < latt_ln ; i++)
	{
		for(int j = 0 ; j < latt_ln ; j++)
		{
			phi_old[0][i][j] = phi[0][i][j];
			phi_old[1][i][j] = phi[1][i][j];
		}
	}
}

long double toda_hybrid_3::hamil()
{
	long double h = 0;
	int ln_m = latt_ln - 1;

	for(int i = 0 ; i < latt_ln ; i++)
	{
		for(int j = 0 ; j < latt_ln ; j++)
		{
			h += 0.5 * sq(pi[0][i][j]) +

			     0.5 * sq(pi[1][i][j]) +

			     0.5 * (sq(phi[0][i][j] - phi[0][(i+1) % latt_ln][j]) +

				        sq(phi[0][i][j] - phi[0][i][(j + 1) % latt_ln])) +
	   
	             0.5 * (sq(phi[1][i][j] - phi[1][(i+1) % latt_ln][j]) +

				        sq(phi[1][i][j] - phi[1][i][(j + 1) % latt_ln])) +

	             c[0] * exp(beta_alpha[0][0] * phi[0][i][j] + beta_alpha[0][1] * phi[1][i][j]) +

                 c[1] * exp(beta_alpha[1][0] * phi[0][i][j] + beta_alpha[1][1] * phi[1][i][j]) +

                 c[2] * exp(beta_alpha[2][0] * phi[0][i][j] + beta_alpha[2][1] * phi[1][i][j]);

		}
	}
	return h;
}

void toda_hybrid_3::pi_hlf_stp()
{
	int ln_m = latt_ln - 1;
	for(int i = 0 ; i < latt_ln ; i++)
	{
		for(int j = 0 ; j < latt_ln ; j++)
		{
			pi[0][i][j] =  pi[0][i][j] + 

			             (((-phi[0][i][j] + phi[0][(i+1) % latt_ln][j]) +

				          (-phi[0][i][j] + phi[0][(i + ln_m) % latt_ln][j]) +

				          (-phi[0][i][j] + phi[0][i][(j + 1) % latt_ln]) +

		                  (-phi[0][i][j] + phi[0][i][(j + ln_m) % latt_ln])) - 

			           	  c[0] * beta_alpha[0][0] * exp(beta_alpha[0][0] * phi[0][i][j] + beta_alpha[0][1] * phi[1][i][j]) -

               			  c[1] * beta_alpha[1][0] * exp(beta_alpha[1][0] * phi[0][i][j] + beta_alpha[1][1] * phi[1][i][j]) -

                		  c[2] * beta_alpha[2][0] * exp(beta_alpha[2][0] * phi[0][i][j] + beta_alpha[2][1] * phi[1][i][j])

			               ) * (tau / 2.0);

			pi[1][i][j] =  pi[1][i][j] + 

			             (((-phi[1][i][j] + phi[1][(i+1) % latt_ln][j]) +

				          (-phi[1][i][j] + phi[1][(i + ln_m) % latt_ln][j]) +

				          (-phi[1][i][j] + phi[1][i][(j + 1) % latt_ln]) +

		                  (-phi[1][i][j] + phi[1][i][(j + ln_m) % latt_ln])) -

			              c[0] * beta_alpha[0][1] * exp(beta_alpha[0][0] * phi[0][i][j] + beta_alpha[0][1] * phi[1][i][j]) -

               			  c[1] * beta_alpha[1][1] * exp(beta_alpha[1][0] * phi[0][i][j] + beta_alpha[1][1] * phi[1][i][j]) -

                		  c[2] * beta_alpha[2][1] * exp(beta_alpha[2][0] * phi[0][i][j] + beta_alpha[2][1] * phi[1][i][j])

			               ) * (tau / 2.0);
		}
	}
}

void toda_hybrid_3::pi_fll_stp()
{
	int ln_m = latt_ln - 1;
	for(int i = 0 ; i < latt_ln ; i++)
	{
		for(int j = 0 ; j < latt_ln ; j++)
		{
			pi[0][i][j] =  pi[0][i][j] + 

			             (((-phi[0][i][j] + phi[0][(i+1) % latt_ln][j]) +

				          (-phi[0][i][j] + phi[0][(i + ln_m) % latt_ln][j]) +

				          (-phi[0][i][j] + phi[0][i][(j + 1) % latt_ln]) +

		                  (-phi[0][i][j] + phi[0][i][(j + ln_m) % latt_ln])) -

			           	  c[0] * beta_alpha[0][0] * exp(beta_alpha[0][0] * phi[0][i][j] + beta_alpha[0][1] * phi[1][i][j]) -

               			  c[1] * beta_alpha[1][0] * exp(beta_alpha[1][0] * phi[0][i][j] + beta_alpha[1][1] * phi[1][i][j]) -

                		  c[2] * beta_alpha[2][0] * exp(beta_alpha[2][0] * phi[0][i][j] + beta_alpha[2][1] * phi[1][i][j])

			               ) * tau;

			pi[1][i][j] =  pi[1][i][j] + 

			             (((-phi[1][i][j] + phi[1][(i+1) % latt_ln][j]) +

				          (-phi[1][i][j] + phi[1][(i + ln_m) % latt_ln][j]) +

				          (-phi[1][i][j] + phi[1][i][(j + 1) % latt_ln]) +

		                  (-phi[1][i][j] + phi[1][i][(j + ln_m) % latt_ln])) -

			              c[0] * beta_alpha[0][1] * exp(beta_alpha[0][0] * phi[0][i][j] + beta_alpha[0][1] * phi[1][i][j]) -

               			  c[1] * beta_alpha[1][1] * exp(beta_alpha[1][0] * phi[0][i][j] + beta_alpha[1][1] * phi[1][i][j]) -

                		  c[2] * beta_alpha[2][1] * exp(beta_alpha[2][0] * phi[0][i][j] + beta_alpha[2][1] * phi[1][i][j])

			               ) * tau;
		}
	}
}

bool toda_hybrid_3::update()
{
	
	old_to_phi();
    pi_init();
    long double eb = hamil();
    pi_hlf_stp();
    phi_fll_stp();
    for(int i = 0 ; i < nmd - 1 ; i++)
    {
        pi_fll_stp();
        phi_fll_stp();
    }
    pi_hlf_stp();
    long double ea = hamil();
    bool flag;
    
    if(ea <= eb || exp(eb - ea) > gsl_rng_uniform(rlx))
    {
        flag = true;
    }
    else
    {
        phi_to_old();
        flag = false;
    }
    return flag;
}


void toda_hybrid_3::phi_fll_stp()
{
	int ln_m = latt_ln - 1;
	for(int i = 0 ; i < latt_ln ; i++)
	{
		for(int j = 0 ; j < latt_ln ; j++)
		{
			phi[0][i][j] =  phi[0][i][j] + pi[0][i][j] * tau;
			phi[1][i][j] =  phi[1][i][j] + pi[1][i][j] * tau;
		}
	}
}

void toda_hybrid_3::smp_mn_fld()
{
	for(int i = 0 ; i < latt_ln ; i++)
	{
		for(int j = 0 ; j < latt_ln ; j++)
		{
			mn_fld[0] += phi[0][i][j];
			mn_fld[1] += phi[1][i][j];
		}
	}
}


void toda_hybrid_3::smp_wall_corr_buff(int N)
{
	vector<long double> fldmnho0(latt_ln , 0.0);
	vector<long double> fldmnve0(latt_ln , 0.0);
	vector<long double> fldmnho1(latt_ln , 0.0);
	vector<long double> fldmnve1(latt_ln , 0.0);

	
	for(int i = 0 ; i < latt_ln ; i++)
	{
		for(int j = 0 ; j < latt_ln ;  j++)
		{
			fldmnve0[i] += phi[0][j][i];
			fldmnve1[i] += phi[1][j][i];
			fldmnho0[i] += phi[0][i][j];
			fldmnho1[i] += phi[1][i][j];
		}
	}

	for(int i = 0 ; i < fit_ln ; i++)
	{
		
		for(int j = 0 ; j < latt_ln ;  j++)
		{
			wall_corr_buff[i][N / num_smps_per_buff] += fldmnve0[j] * fldmnve0[(j + i) % latt_ln]
			                                          + fldmnho0[j] * fldmnho0[(j + i) % latt_ln];

			wall_corr_buff[i + fit_ln][N / num_smps_per_buff] += fldmnve1[j] * fldmnve1[(j + i) % latt_ln]
								                               + fldmnho1[j] * fldmnho1[(j + i) % latt_ln];

			wall_corr_buff[i+ 2 * fit_ln][N / num_smps_per_buff] += 0.5 * (fldmnve0[j] * fldmnve1[(j + i) % latt_ln]
			                                                             + fldmnho0[j] * fldmnho1[(j + i) % latt_ln]
			                                                             + fldmnve1[j] * fldmnve0[(j + i) % latt_ln]
			                                                             + fldmnho1[j] * fldmnho0[(j + i) % latt_ln]);

		}
	}
	
}

void toda_hybrid_3::calc_wall_corr_frm_buff()
{
	for(int i = 0 ; i < 3 * fit_ln ; i++)
	{
		for(int j = 0 ; j < num_buffs ; j++)
		{
			wall_corr[i] += wall_corr_buff[i][j];
		}
	}
}


void toda_hybrid_3::crt_bnd_dat(int  num_buffs_per_bin_log2)
{
	this -> num_buffs_per_bin_log2 = num_buffs_per_bin_log2;
	num_buffs_per_bin = 1 << num_buffs_per_bin_log2;
	num_bins = num_buffs / num_buffs_per_bin;
	bnd_dat.resize(3 * fit_ln , vector<long double>(num_bins, 0.0));

	for(int i = 0 ; i < 3 * fit_ln ; i++)
	{
		for(int j = 0 ; j < num_bins ; j++)
		{
			for(int k = 0 ; k < num_buffs_per_bin ; k++)
			{
				bnd_dat[i][j] += wall_corr_buff[i][j * num_buffs_per_bin + k];
			}

			bnd_dat[i][j] /= (num_buffs_per_bin * num_smps_per_buff);
		}
	}
}

void toda_hybrid_3::calc_wall_corr_frm_bnd_dat(vector<long double> &test_corr)
{
	
	long long div = 2 * latt_ln * latt_ln * latt_ln * ((long long)num_bins);
   
	for(int i = 0 ; i < 3 * fit_ln ; i++)
	{
		for(int j = 0 ; j < num_bins ; j++)
		{
			test_corr[i] += bnd_dat[i][j];
		}
		test_corr[i] /= div;
		
	}
}


void toda_hybrid_3::crt_btstrp_wall_corr(double *w_corr_fit)
{
	
	for(int i = 0 ; i < 3 * fit_ln ; i++)
	{
		w_corr_fit[i] = 0.0;
	}

    int ind;

	for(int j = 0 ; j < num_bins ; j++)
	{
		ind = gsl_rng_uniform(rlx) * num_bins;
		for(int i = 0 ; i < 3 * fit_ln ; i++)
		{
			w_corr_fit[i] += bnd_dat[i][ind];
		}
	}

	long long div = 2 * latt_ln * latt_ln * latt_ln * ((long long)num_bins);

	for(int i = 0 ; i < 3 * fit_ln ; i++)
	{
		w_corr_fit[i] /= div;
	}
}


void toda_hybrid_3::calc_err(double max_err[])
{
	double  temp[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    for(int i = 0 ; i < 9 ; i++)
    {
    	max_err[i] = 0.0;
    }

	for(int i = 2 ; i < (iters_log2 - num_smps_per_buff_log_2) ; i++)
    {
        calc_err_at_bn_dpth(i,5000, temp);
        for(int j = 0 ; j < 9 ; j++)
        {
        	if(temp[j] > max_err[j])
	        {
	        	max_err[j] = temp[j];
	        } 
        }
    }
    
}

void toda_hybrid_3::calc_err_at_bn_dpth(int bin_dpth, int num_bstp_smps, double err[])
{
	int Nfit = 3 * fit_ln;
    double y[Nfit];
    double x_init[8] = {2.0, 2.0, 2.0, ((double)mn_fld[0]), 2.0, 2.0, 2.0, ((double)mn_fld[1])};
    double res[9];
	crt_bnd_dat(bin_dpth);

	double sum_sum_sq[2][9];
	for (int i = 0 ; i < 9 ; i++)
	{
		sum_sum_sq[0][i] = 0.0;
		sum_sum_sq[1][i] = 0.0;
	}

	for(int i = 0 ; i < 100 ; i++)
	{
		crt_btstrp_wall_corr(y);
        mass_ratio_fit(Nfit, y, x_init, res);
	}

    int Nbt = num_bstp_smps;
    for(int i = 0 ; i < Nbt ; i++)
    {
        crt_btstrp_wall_corr(y);
        
        mass_ratio_fit(Nfit, y, x_init, res);
        for(int j = 0 ; j < 9 ; j++)
        {
        	sum_sum_sq[0][j] += res[j];
        	sum_sum_sq[1][j] += res[j] * res[j];
        }
        
    }
    for (int i = 0 ; i < 9 ; i++)
	{
		err[i] = (sum_sum_sq[1][i]/Nbt) - (sum_sum_sq[0][i]/Nbt) * (sum_sum_sq[0][i]/Nbt);
	}
    
}




void toda_hybrid_3::prnt_mn_fld()
{
	cout << "mean phi_0 field = " << mn_fld[0]/(latt_ln * latt_ln * iters) << endl;
	cout << "mean phi_1 field = " << mn_fld[1]/(latt_ln * latt_ln * iters) << endl;
}

void toda_hybrid_3::prnt_results()
{
	string file_out = type;

    file_out +=  "_lattln" + to_string(latt_ln) + "_beta" + to_string(beta) +
	                  "_mass" + to_string(m) + "_tau" + to_string(tau)  + "_nmd"
	                  + to_string(nmd) + "_iters" + to_string(iters) + "_fitln" + to_string(fit_ln) + ".dat";

	ofstream write_output(file_out);
	assert(write_output.is_open());

    mn_fld[0] /= (latt_ln * latt_ln * ((long long)iters));
    mn_fld[1] /= (latt_ln * latt_ln * ((long long)iters));

	write_output << "mean phi_0 field = " << mn_fld[0] << endl;
	write_output << "mean phi_1 field = " << mn_fld[1] << endl;

	write_output << "acceptance_rate = " << ((double)iters)/attempts << endl;
    write_output << "corrected version" << endl << endl;

	long long  div = 2 * latt_ln * latt_ln * latt_ln * ((long long)iters);
    calc_wall_corr_frm_buff();
	for(int i = 0 ; i < 3 * fit_ln ; i++)
	{
		wall_corr[i] /= div;
	}

	write_output << "corr = [";
	for(int i = 0 ; i < 3 * fit_ln - 1 ; i++)
	{
		write_output << wall_corr[i] << ", "; 
	}
    write_output << wall_corr[3 * fit_ln - 1] << "];" << endl << endl;
    
    size_t Nfit = 3 * fit_ln;
    double y[Nfit];
    for(int i = 0 ; i < Nfit ; i++)
    {
        y[i] = (double)wall_corr[i];
    }
    double x_init[8] = {2.0,2.0,2.0, (double)mn_fld[0],2.0,2.0,2.0, (double)mn_fld[1]};
    
    double results[9];
    mass_ratio_fit(Nfit, y, x_init, results);

    double err[9];
    calc_err(err);
    
    for (int i = 0 ; i < 9 ; i++)
    {
        write_output << results[i] << "  " << err[i] << endl;
    }
    

    write_output.close();
}


toda_hybrid_3::~toda_hybrid_3()
{
	gsl_rng_free(rlx);
	gsl_rng_free(r48);
}