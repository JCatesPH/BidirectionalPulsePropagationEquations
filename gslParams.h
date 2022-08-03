#pragma once

using namespace std;

#include <fftw3.h>
#include <complex> 
#include <cmath>
#include <map>
#include "gsl/gsl_math.h"

#include "Materials.h"

#define FFTW_WISDOM_TYPE FFTW_PATIENT

class ODEParams { // Class for ODE params
	private:
		/* map<string, gsl_odeiv2_step_type> odeStepMap {{"RK4", *gsl_odeiv2_step_rk4}, 
			{"RK45_Fehlberg", *gsl_odeiv2_step_rkf45},
			{"RK45_CashKarp", *gsl_odeiv2_step_rkck},
			{"RK45_CashKarp", *gsl_odeiv2_step_rkck},
			{"RK89_PrinceDormand", *gsl_odeiv2_step_rk8pd},
			{"MSADAMS", *gsl_odeiv2_step_msadams},
		}; */
    public:
		int numT = -666;
		int numX = -666;
		int numOmeg = -666;
		int numOmX = -666;
		int doPlasmaCalc = 0;
		double chi_2, chi_3, mpi_sigmaK, mpi_k, ionE, sigmaBremsstrahlung, recombTime;
		double *omega, *kx, *rho, *y;
		complex<double> *k, *ee_p, *ee_m, *p_nl, *jhat, *j_e;
		fftw_plan p_ffft, j_ffft, ep_b, em_b, ep_f, em_f;
		int printAtomicProfile = 0;
		//gsl_odeiv2_step_type stepType;

		ODEParams() { // Default constructor (empty)
			omega = nullptr;
			kx = nullptr;
			rho = nullptr;
			y = nullptr;
			k = nullptr;
			ee_p = nullptr;
			ee_m = nullptr;
			p_nl = nullptr;
			jhat = nullptr;
			j_e = nullptr;
			p_ffft = nullptr;
			j_ffft = nullptr;
			ep_b = nullptr;
			em_b = nullptr;
			ep_f = nullptr;
			em_f = nullptr;
		}

		
		ODEParams(int nT, double *omg); // Constructor definition used most often

		ODEParams(int nT, int nX, double *omg, double* kperp);

		~ODEParams() { // Destructor definition
			// Free arrays in heap
			free(rho);
			free(y);
			free(ee_p);
			free(ee_m);
			free(p_nl);
			free(jhat);
			free(j_e);

			// Destroy fftw plans
			fftw_destroy_plan(p_ffft);
			fftw_destroy_plan(j_ffft);
			fftw_destroy_plan(ep_b);
			fftw_destroy_plan(em_b);
			fftw_destroy_plan(ep_f);
			fftw_destroy_plan(em_f);
		}

		// Copy constructor
		ODEParams(const ODEParams &p1);

		// Function for filling params with structure properties
		void fillParams(Material myMat);

		void initializeY(complex<double> *sL);

		/* void setStepType(string inStr) { stepType = odeStepMap[inStr];}
		gsl_odeiv2_step_type getStepType() { return stepType;} */
};




/*===========================================================================*/

class RootParams {
	private:
		int m_itnum = 1;
		int m_output = 0;
		int m_intCondition = 0; 
		int m_nRoot = -666;
		double m_Gnorm;
		ODEParams *m_odeObj = nullptr;

		complex<double> *m_sourceLeft, *m_linearRefl;
		double m_lagrange = 0.0;
		double *derivFactor;
    public:
		complex<double> *integral;
		RootParams() {} // Default constructor

      	RootParams(ODEParams *odeObj, int sizeRoot) { // Constructor Definition
			m_nRoot = sizeRoot;
			m_odeObj = odeObj;
			//m_odeObj = ODEParams(nT, omg);
			derivFactor = (double*)malloc(sizeof(double) * 2 * sizeRoot);

			for (int j = 0; j < 2 * sizeRoot; j++) {
				derivFactor[j] = GSL_SQRT_DBL_EPSILON;
			}
		} 

		/* RootParams(int nT, double *omg, int sizeRoot) { // Constructor Definition
			m_nRoot = sizeRoot;
			//m_odeObj = odeObj;
			m_odeObj = new ODEParams;
			*m_odeObj = ODEParams(nT, omg);
		} */

		RootParams(const RootParams &root_in) { // Copy constructor
			m_itnum = root_in.getItNum();
			m_output = root_in.getOutParam();
			m_intCondition = root_in.getIntCond();
			m_nRoot = root_in.getSizeRoot();

			printf("     Copying ODEParams in copy of RootParams\n");
			//ODEParams odeObjNew = *(root_in.m_odeObj);
			m_odeObj = new ODEParams;
			//m_odeObj->cloneFrom(root_in.getODEparams());
			*m_odeObj = *(root_in.getODEparams());

			m_sourceLeft = root_in.getSourceLeft();
			m_linearRefl = root_in.getLinSol();
		}

		~RootParams() { // Destructor definition
			free(derivFactor);
			//m_odeObj->~ODEParams();
		}

		int getItNum() const { return m_itnum; }
		int getOutParam() const { return m_output; }
		int getIntCond() const { return m_intCondition; }
		int getSizeRoot() const { return m_nRoot; }
		ODEParams *getODEparams() const { return m_odeObj; }
		double getGnorm() const { return m_Gnorm; }
		complex<double> * getSourceLeft() const { return m_sourceLeft; }
		complex<double> * getLinSol() const { return m_linearRefl; }
		double getLagrangeMult() const { return m_lagrange; }

		void setItNum(int aItNum) { m_itnum = aItNum; }
		void setOutParam(int aOutParam) { m_output = aOutParam; }
		void setIntCond(int aIntCond) { m_intCondition = aIntCond; }
		void setODEparams(ODEParams *aODEobj) { m_odeObj = aODEobj; }
		void setGnorm(double fout) { m_Gnorm = fout; }
		void setSourceLeft(complex<double> *sIn) { m_sourceLeft = sIn; }
		void setLinSol(complex<double> *aIn) { m_linearRefl = aIn; }
		void setLagrangeMult(double mu) { m_lagrange = mu; }
		double *getDerivFactors() const { return derivFactor; }
};
