#pragma once

using namespace std;

#include <fftw3.h>
#include <complex> 
#include <cmath>

typedef struct {
    int numT, numOmeg, doPlasmaCalc;
	double chi_2, chi_3, mpi_sigmaK, mpi_k, ionE;
	double *omega, *kx, *rho, *y;
	complex<double> *k, *ee_p, *ee_m, *nl_k, *nl_p, *j_e;
	fftw_plan nk_f, ep_b, em_b, np_f, ep_f, em_f;
}odeparam_type;

typedef struct {
	int itnum, output, intCondition, nRoot;
	odeparam_type *odestruct;
}rootparam_type;


odeparam_type* fill_odeparams(
    int nT,
    double chi_2, 
    double chi_3, 
    double* omg, 
    double* kx,
    double* nE,
    complex<double> *j_e,
    complex<double>* k, 
    complex<double> *eForward,
    complex<double> *eBackward,
    complex<double> *nonlinearP,
    complex<double> *nlJhat,
    double *yp);

odeparam_type* copy_odeparams(odeparam_type *p_in);
rootparam_type *copy_rootparams(rootparam_type *rp_in);
void free_odeparams(odeparam_type * p_in);
void free_rootparams(rootparam_type *rp_in);