#include "mass_curve_fit.hpp"

#define PI 3.14159265359

struct data 
{
	size_t n;
	double * t;
	double * y;
};

int expb_f (const gsl_vector * x , void *data,
                                gsl_vector * f)
{
    size_t n = ((struct data *)data)->n;
    double *t = ((struct data *)data)->t;
    double *y = ((struct data *)data)->y;

    double a = gsl_vector_get (x, 0);
    double theta_1 = gsl_vector_get (x, 1);
    double m_1 = gsl_vector_get (x, 2);
    double S_1 = gsl_vector_get (x, 3);
    double b = gsl_vector_get (x, 4);
    double theta_2 = gsl_vector_get (x, 5);
    double m_2 = gsl_vector_get (x, 6);
    double S_2 = gsl_vector_get (x, 7);

    int fit_ln = n/3;
   /*
    double t[n];
    int fit_ln = n/3;
    for (int i = 0; i < n; i++)
    {
        t[i] = (i%fit_ln)*1.0;
    }*/

    for (size_t i = 0; i < n; i++)
    {  
		double Yi =

		(i < fit_ln) * (a * (cos(theta_1) * cos(theta_1)) * exp(-m_1 * t[i]) +
		                b * (sin(theta_2) * sin(theta_2)) * exp(-m_2 * t[i]) + S_1 * S_1) +

		(i >= fit_ln && i < 2 * fit_ln) * (a * (sin(theta_1) * sin(theta_1)) * exp(-m_1 * t[i]) +
		                                   b * (cos(theta_2) * cos(theta_2)) * exp(-m_2 * t[i]) + S_2 * S_2) +

		(i >= 2 * fit_ln && i < 3 * fit_ln) * (-a * (cos(theta_1) * sin(theta_1)) * exp(-m_1 * t[i])+
		                                        b * (sin(theta_2) * cos(theta_2)) * exp(-m_2 * t[i]) + S_1 * S_2);
		gsl_vector_set (f, i, Yi - y[i]);
    }

    return GSL_SUCCESS;
}

/*************************************************************************************/

int expb_df (const gsl_vector * x, void *params, gsl_matrix * J)
{
    size_t n = ((struct data *)params)->n;
    double *t = ((struct data *)params)->t;

    double a = gsl_vector_get (x, 0);
    double theta_1 = gsl_vector_get (x, 1);
    double m_1 = gsl_vector_get (x, 2);
    double S_1 = gsl_vector_get (x, 3);
    double b = gsl_vector_get (x, 4);
    double theta_2 = gsl_vector_get (x, 5);
    double m_2 = gsl_vector_get (x, 6);
    double S_2 = gsl_vector_get (x, 7);
    
    int fit_ln = n/3;
    /*double t[n];
    int fit_ln = n/3;

    for (int i = 0; i < n; i++)
    {
        t[i] = (i % fit_ln)*1.0;
    }*/

    double J_0, J_1, J_2, J_3, J_4, J_5, J_6, J_7;
    
	for (size_t i = 0; i < n; i++)
	{
		J_0 = (i < fit_ln) * ((cos(theta_1) * cos(theta_1)) * exp(-m_1 * t[i])) +

		      (i >= fit_ln && i < 2 * fit_ln) * ((sin(theta_1) * sin(theta_1)) * exp(-m_1 * t[i])) +

		      (i >= 2 * fit_ln && i < 3 * fit_ln) * (-(cos(theta_1) * sin(theta_1)) * exp(-m_1 * t[i]));

		J_1 = (i < fit_ln) * ( -2 * a * (sin(theta_1) * cos(theta_1)) * exp( -m_1 * t[i])) +

		      (i >= fit_ln && i < 2 * fit_ln) * (2 * a * (cos(theta_1) * sin(theta_1)) * exp(-m_1 * t[i])) +

		      (i >= 2 * fit_ln && i < 3 * fit_ln) * (-a * ((cos(theta_1) * cos(theta_1) - 
		                                                    sin(theta_1) * sin(theta_1)) * exp(-m_1 * t[i])));

		J_2 = (i < fit_ln) * (a * (cos(theta_1) * cos(theta_1)) * -t[i] * exp(-m_1 * t[i])) +

		      (i >= fit_ln && i < 2 * fit_ln) * (a * (sin(theta_1) * sin(theta_1)) * -t[i] * exp(-m_1 * t[i])) +

		      (i >= 2 * fit_ln && i < 3 * fit_ln) * (a * (cos(theta_1) * sin(theta_1)) *  t[i] * exp(-m_1 * t[i]));

		J_3 = (i < fit_ln) * (2 * S_1) +

		      (i >= fit_ln && i < 2 * fit_ln) * (0) +

		      (i >= 2 * fit_ln && i < 3 * fit_ln) * (S_2);

		J_4 = (i < fit_ln) * ((sin(theta_2) * sin(theta_2)) * exp(-m_2 * t[i])) +

		      (i >= fit_ln && i < 2 * fit_ln) * ((cos(theta_2) * cos(theta_2)) * exp(-m_2 * t[i])) +

		      (i >= 2 * fit_ln && i < 3 * fit_ln) * ((sin(theta_2) * cos(theta_2)) * exp(-m_2 * t[i]));

		J_5 = (i < fit_ln) * (2 * b * (cos(theta_2) * sin(theta_2)) * exp(-m_2 * t[i])) +

		      (i >= fit_ln && i < 2 * fit_ln) * ( -2 * b * (sin(theta_2) * cos(theta_2)) * exp(-m_2 * t[i])) +

		      (i >= 2 * fit_ln && i < 3 * fit_ln) * ( b * (cos(theta_2) * cos(theta_2)
		                            - sin(theta_2) * sin(theta_2)) * exp(-m_2 * t[i]));

		J_6 = (i < fit_ln) * b * (sin(theta_2) * sin(theta_2)) * -t[i] * exp(-m_2 * t[i]) +

		      (i >= fit_ln && i < 2 * fit_ln) * (b * (cos(theta_2) * cos(theta_2)) * -t[i] * exp(-m_2 * t[i])) +

		      (i >= 2 * fit_ln && i < 3 * fit_ln) * (b * (sin(theta_2) * cos(theta_2)) * -t[i] * exp(-m_2 * t[i]));

		J_7 = (i < fit_ln) * (0) + 

		      (i >= fit_ln && i < 2 * fit_ln) * (2 * S_2) +
		      
              (i >= 2 * fit_ln && i < 3 * fit_ln) * (S_1);

		gsl_matrix_set(J, i, 0, J_0);
		gsl_matrix_set(J, i, 1, J_1);
		gsl_matrix_set(J, i, 2, J_2);
		gsl_matrix_set(J, i, 3, J_3);
		gsl_matrix_set(J, i, 4, J_4);
		gsl_matrix_set(J, i, 5, J_5);
		gsl_matrix_set(J, i, 6, J_6);
		gsl_matrix_set(J, i, 7, J_7);
	}	
    return GSL_SUCCESS;
}

void
callback(const size_t iter, void *params,
         const gsl_multifit_nlinear_workspace *w)
{
  
}


/*************************************************************************************/

void mass_ratio_fit(size_t n, double y[], double x_init[], double res[])
{
	const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
  	gsl_multifit_nlinear_workspace *w;
  	gsl_multifit_nlinear_fdf fdf;
  	gsl_multifit_nlinear_parameters fdf_params =
    gsl_multifit_nlinear_default_parameters();

    const size_t p = 8;
    gsl_vector *f;
  	gsl_matrix *J;
  	gsl_matrix *covar = gsl_matrix_alloc (p, p);
  	double t[n], weights[n];
 
    int fit_ln = n/3;

    for (int i = 0; i < n; i++)
    {
        weights[i] = 1.0;
        t[i] = (i % fit_ln)*1.0;
    }


  	struct data d = { n, t, y };

  	gsl_vector_view x = gsl_vector_view_array (x_init, p);
  	gsl_vector_view wts = gsl_vector_view_array(weights, n);

    double chisq, chisq0;
	int status, info;
	size_t i;

	const double xtol = 1e-8;
	const double gtol = 1e-8;
	const double ftol = 0.0;

    fdf.f = expb_f;
  	fdf.df = expb_df;   /* set to NULL for finite-difference Jacobian */
  	fdf.fvv = NULL;     /* not using geodesic acceleration */
  	fdf.n = n;
  	fdf.p = p;
  	fdf.params = &d;

    /* allocate workspace with default parameters */
  	w = gsl_multifit_nlinear_alloc (T, &fdf_params, n, p);

  	/* initialize solver with starting point and weights */
 	gsl_multifit_nlinear_winit (&x.vector, &wts.vector, &fdf, w);



    status = gsl_multifit_nlinear_driver(100, xtol, gtol, ftol,
                                       callback, NULL, &info, w);

	size_t iter = 0;

	
    
    for(int i = 0 ; i < 8 ; i++)
    {
    	res[i] = gsl_vector_get(w->x, i);
    }
	double ratio_1 = (gsl_vector_get(w->x, 2))/(gsl_vector_get(w->x, 6));
	double ratio_2 = (gsl_vector_get(w->x, 6))/(gsl_vector_get(w->x, 2));
            
    double mass_ratio;

	if (ratio_1 > ratio_2) 
	{
		res[8] = ratio_1;
	}
	else
	{
		res[8] = ratio_2;
	}

	res[1] = fmod(res[1], 2 * PI);
	res[5] = fmod(res[5], 2 * PI);
    
    

	gsl_multifit_nlinear_free (w);
  gsl_matrix_free (covar);
            
	  

}
