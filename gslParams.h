#pragma once

using namespace std;

#include <fftw3.h>
#include <complex> 
#include <cmath>

#include "Materials.h"

#define FFTW_WISDOM_TYPE FFTW_PATIENT

class ODEParams { // Class for ODE params
    public:
		int numT, numOmeg, doPlasmaCalc;
		double chi_2, chi_3, mpi_sigmaK, mpi_k, ionE;
		double *omega, *kx, *rho, *y;
		complex<double> *k, *ee_p, *ee_m, *nl_k, *nl_p, *j_e;
		fftw_plan nk_f, ep_b, em_b, np_f, ep_f, em_f;

		//ODEParams() {} // Default constructor
		
		ODEParams(int nT, double *omg) { // Constructor definition
			// Takes parameters to set values that are global.
			numT = nT;
			numOmeg = numT / 2 + 1;
			omega = omg;

			// Allocates necessary double arrays
			rho = (double*)malloc(sizeof(double)*numT);
        	y = (double*)malloc(sizeof(double)*(4*numOmeg));

			// Allocate complex arrays
			ee_p = (complex<double>*)malloc(sizeof(complex<double>)*numT);
			ee_m = (complex<double>*)malloc(sizeof(complex<double>)*numT);
			nl_k = (complex<double>*)malloc(sizeof(complex<double>)*numT);
			nl_p = (complex<double>*)malloc(sizeof(complex<double>)*numT);
			j_e = (complex<double>*)malloc(sizeof(complex<double>)*numT);

			// Allocate fftw plans
			nk_f = fftw_plan_dft_1d(numT, reinterpret_cast<fftw_complex*>(&nl_k[0]), reinterpret_cast<fftw_complex*>(&nl_k[0]), FFTW_FORWARD, FFTW_WISDOM_TYPE );
			np_f = fftw_plan_dft_1d(numT, reinterpret_cast<fftw_complex*>(&j_e[0]), reinterpret_cast<fftw_complex*>(&nl_p[0]), FFTW_FORWARD, FFTW_WISDOM_TYPE );
			ep_b = fftw_plan_dft_1d(numT, reinterpret_cast<fftw_complex*>(&ee_p[0]), reinterpret_cast<fftw_complex*>(&ee_p[0]), FFTW_BACKWARD, FFTW_WISDOM_TYPE );
			em_b = fftw_plan_dft_1d(numT, reinterpret_cast<fftw_complex*>(&ee_m[0]), reinterpret_cast<fftw_complex*>(&ee_m[0]), FFTW_BACKWARD, FFTW_WISDOM_TYPE );
			ep_f = fftw_plan_dft_1d(numT, reinterpret_cast<fftw_complex*>(&ee_p[0]), reinterpret_cast<fftw_complex*>(&ee_p[0]), FFTW_FORWARD, FFTW_WISDOM_TYPE );
			em_f = fftw_plan_dft_1d(numT, reinterpret_cast<fftw_complex*>(&ee_m[0]), reinterpret_cast<fftw_complex*>(&ee_m[0]), FFTW_FORWARD, FFTW_WISDOM_TYPE );
		}

		~ODEParams() { // Destructor definition
			// Free arrays in heap
			free(rho);
			free(y);
			free(ee_p);
			free(ee_m);
			free(nl_k);
			free(nl_p);
			free(j_e);

			// Destroy fftw plans
			fftw_destroy_plan(nk_f);
			fftw_destroy_plan(ep_b);
			fftw_destroy_plan(em_b);
			fftw_destroy_plan(np_f);
			fftw_destroy_plan(ep_f);
			fftw_destroy_plan(em_f);
		}

		// Copy constructor
		ODEParams(const ODEParams &p1);

		// Function for filling params with structure properties
		void fillParams(Material myMat);
};




/*===========================================================================*/

class RootParams {
	private:
		int m_itnum = 1;
		int m_output = 0;
		int m_intCondition = 0; 
		int m_nRoot;
		ODEParams *m_odeObj;
    public:
		complex<double> *integral;
		RootParams() {} // Default constructor

      	RootParams(ODEParams *odeObj, int sizeRoot) { // Constructor Definition
				m_nRoot = sizeRoot;
				m_odeObj = odeObj;
				//m_odeObj = ODEParams(nT, omg);
			}

		RootParams(const RootParams &p1) { // Copy constructor
			m_itnum = p1.m_itnum;
			m_output = p1.m_output;
			m_intCondition = p1.m_intCondition;
			m_nRoot = p1.m_nRoot;
			printf("     Copying ODEParams in copy of RootParams\n");
			*m_odeObj = *p1.m_odeObj;
		}

		int getItNum() { return m_itnum; }
		int getOutParam() { return m_output; }
		int getIntCond() { return m_intCondition; }
		int getSizeRoot() { return m_nRoot; }
		ODEParams *getODEparams() { return m_odeObj; }

		void setItNum(int aItNum) { m_itnum = aItNum; }
		void setOutParam(int aOutParam) { m_output = aOutParam; }
		void setIntCond(int aIntCond) { m_intCondition = aIntCond; }
		void setODEparams(ODEParams *aODEobj) { m_odeObj = aODEobj; }
};