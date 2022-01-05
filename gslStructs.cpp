#include "gslStructs.h"

#define FFTW_WISDOM_TYPE FFTW_PATIENT

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
    double *yp) {

	odeparam_type* r = (odeparam_type*)malloc(sizeof(odeparam_type));
	if (r != NULL)
	{   
        r->numT = nT;
		r->chi_2 = chi_2;
		r->chi_3 = chi_3;
		r->omega = omg;
		r->kx = kx;
        r->k = k;
		/* r->rho = (double*)malloc(sizeof(double)*nT);
		r->j_e = (complex<double>*)malloc(sizeof(complex<double>)*nT);
		r->ee_p = (complex<double>*)malloc(sizeof(complex<double>)*nT);
		r->ee_m = (complex<double>*)malloc(sizeof(complex<double>)*nT);
		r->nl_k = (complex<double>*)malloc(sizeof(complex<double>)*nT);
		r->nl_p = (complex<double>*)malloc(sizeof(complex<double>)*nT); */
        r->rho = nE;
        r->j_e = j_e;
        r->ee_p = eForward;
        r->ee_m = eBackward;
        r->nl_k = nonlinearP;
        r->nl_p = nlJhat;
		r->nk_f = fftw_plan_dft_1d(r->numT, reinterpret_cast<fftw_complex*>(&r->nl_k[0]), reinterpret_cast<fftw_complex*>(&r->nl_k[0]), FFTW_FORWARD, FFTW_WISDOM_TYPE );
		r->np_f = fftw_plan_dft_1d(r->numT, reinterpret_cast<fftw_complex*>(&r->j_e[0]), reinterpret_cast<fftw_complex*>(&r->nl_p[0]), FFTW_FORWARD, FFTW_WISDOM_TYPE );
		r->ep_b = fftw_plan_dft_1d(r->numT, reinterpret_cast<fftw_complex*>(&r->ee_p[0]), reinterpret_cast<fftw_complex*>(&r->ee_p[0]), FFTW_BACKWARD, FFTW_WISDOM_TYPE );
		r->em_b = fftw_plan_dft_1d(r->numT, reinterpret_cast<fftw_complex*>(&r->ee_m[0]), reinterpret_cast<fftw_complex*>(&r->ee_m[0]), FFTW_BACKWARD, FFTW_WISDOM_TYPE );
		r->ep_f = fftw_plan_dft_1d(r->numT, reinterpret_cast<fftw_complex*>(&r->ee_p[0]), reinterpret_cast<fftw_complex*>(&r->ee_p[0]), FFTW_FORWARD, FFTW_WISDOM_TYPE );
		r->em_f = fftw_plan_dft_1d(r->numT, reinterpret_cast<fftw_complex*>(&r->ee_m[0]), reinterpret_cast<fftw_complex*>(&r->ee_m[0]), FFTW_FORWARD, FFTW_WISDOM_TYPE );
		r->y = yp;
	}
	else {
		printf("Could not maloc space in fill_params()\n");
		exit(-1);
	}
	return r;
}

odeparam_type* copy_odeparams(odeparam_type *p_in) {
    /* Returns a copy of the input odeparam_type */

	odeparam_type* p_out = (odeparam_type*)malloc(sizeof(odeparam_type));
	if (p_out != NULL)
	{
        // Copy the int and double values
        p_out->numT = p_in->numT;
		p_out->chi_2 = p_in->chi_2;
		p_out->chi_3 = p_in->chi_3;
        p_out->mpi_k = p_in->mpi_k;
        p_out->mpi_sigmaK = p_in->mpi_sigmaK;
        p_out->ionE = p_in->ionE;

        // Copy pointers if data is only read.
		p_out->omega = p_in->omega;
		p_out->kx = p_in->kx;
        p_out->k = p_in->k;

        // Allocate double arrays
        int numOmeg = p_out->numT / 2 + 1;
        p_out->rho = (double*)malloc(sizeof(double)*p_out->numT);
        p_out->y = (double*)malloc(sizeof(double)*(4*numOmeg));

        // Allocate complex arrays
        p_out->ee_p = (complex<double>*)malloc(sizeof(complex<double>)*p_out->numT);
        p_out->ee_m = (complex<double>*)malloc(sizeof(complex<double>)*p_out->numT);
        p_out->nl_k = (complex<double>*)malloc(sizeof(complex<double>)*p_out->numT);
        p_out->nl_p = (complex<double>*)malloc(sizeof(complex<double>)*p_out->numT);
        p_out->j_e = (complex<double>*)malloc(sizeof(complex<double>)*p_out->numT);

        for (int i = 0; i < p_out->numT; i++) {
            p_out->rho[i] = p_in->rho[i];
            p_out->y[i] = p_in->y[i];
            p_out->y[i + numOmeg] = p_in->y[i + numOmeg];
            p_out->y[i + 2*numOmeg] = p_in->y[i + 2*numOmeg];
            p_out->y[i + 3*numOmeg] = p_in->y[i + 3*numOmeg];

            p_out->ee_p[i] = p_in->ee_p[i];
		    p_out->ee_m[i] = p_in->ee_m[i];
		    p_out->nl_k[i] = p_in->nl_k[i];
		    p_out->nl_p[i] = p_in->nl_p[i];
		    p_out->j_e[i] = p_in->j_e[i];
        }

        // Allocate fftw_plans
		p_out->nk_f = fftw_plan_dft_1d(p_out->numT, reinterpret_cast<fftw_complex*>(&p_out->nl_k[0]), reinterpret_cast<fftw_complex*>(&p_out->nl_k[0]), FFTW_FORWARD, FFTW_WISDOM_TYPE );
		p_out->np_f = fftw_plan_dft_1d(p_out->numT, reinterpret_cast<fftw_complex*>(&p_out->j_e[0]), reinterpret_cast<fftw_complex*>(&p_out->nl_p[0]), FFTW_FORWARD, FFTW_WISDOM_TYPE );
		p_out->ep_b = fftw_plan_dft_1d(p_out->numT, reinterpret_cast<fftw_complex*>(&p_out->ee_p[0]), reinterpret_cast<fftw_complex*>(&p_out->ee_p[0]), FFTW_BACKWARD, FFTW_WISDOM_TYPE );
		p_out->em_b = fftw_plan_dft_1d(p_out->numT, reinterpret_cast<fftw_complex*>(&p_out->ee_m[0]), reinterpret_cast<fftw_complex*>(&p_out->ee_m[0]), FFTW_BACKWARD, FFTW_WISDOM_TYPE );
		p_out->ep_f = fftw_plan_dft_1d(p_out->numT, reinterpret_cast<fftw_complex*>(&p_out->ee_p[0]), reinterpret_cast<fftw_complex*>(&p_out->ee_p[0]), FFTW_FORWARD, FFTW_WISDOM_TYPE );
		p_out->em_f = fftw_plan_dft_1d(p_out->numT, reinterpret_cast<fftw_complex*>(&p_out->ee_m[0]), reinterpret_cast<fftw_complex*>(&p_out->ee_m[0]), FFTW_FORWARD, FFTW_WISDOM_TYPE );
	}
	else {
		printf("Could not maloc space in copy_odeparams(..)\n");
		exit(-1);
	}
	return p_out;
}


rootparam_type *copy_rootparams(rootparam_type *rp_in) {
    rootparam_type* rp_out = (rootparam_type*)malloc(sizeof(rootparam_type));
    if (rp_out != NULL)
	{
        rp_out->itnum = rp_in->itnum;
        rp_out->output = rp_in->output;
        rp_out->nRoot = rp_in->nRoot;
        rp_out->intCondition = rp_in->intCondition;
        rp_out->odestruct = copy_odeparams(rp_in->odestruct);
    }
    return rp_out;
}

void free_odeparams(odeparam_type * p_in) {
    if (p_in == NULL) {
        return;
    }
    else {
        // Free all arrays
        free(p_in->rho);
        free(p_in->j_e);
        free(p_in->ee_p);
        free(p_in->ee_m);
        free(p_in->nl_k);
        free(p_in->nl_p);
        free(p_in->y);
        

        // Free fftw_plans
        fftw_destroy_plan(p_in->nk_f);
        fftw_destroy_plan(p_in->ep_b);
        fftw_destroy_plan(p_in->em_b);
        fftw_destroy_plan(p_in->np_f);
        fftw_destroy_plan(p_in->ep_f);
        fftw_destroy_plan(p_in->em_f);

        free(p_in);
        return;
    }
}

void free_rootparams(rootparam_type *rp_in) {
    if (rp_in == NULL) {
        return;
    }
    else {
        free_odeparams(rp_in->odestruct);
        free(rp_in);
        return;
    }
}