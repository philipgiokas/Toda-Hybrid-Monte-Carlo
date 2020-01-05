

#ifndef HOME_MASS_FIT_DEF_3
#define HOME_MASS_FIT_DEF_3

#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cassert>
#include <string>

/*************************************************************************************/

int expb_f (const gsl_vector * x , void *data, gsl_vector * f);

int expb_df (const gsl_vector * x, void *params, gsl_matrix * J);

void callback(const size_t iter, void *params,
         const gsl_multifit_nlinear_workspace *w);

void mass_ratio_fit(size_t n, double y[], double x_init[], double res[]);

#endif