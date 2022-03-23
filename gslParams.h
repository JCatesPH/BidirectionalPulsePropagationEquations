#pragma once

using namespace std;

#include <fftw3.h>
#include <complex> 
#include <cmath>

#include "Materials.h"

#define FFTW_WISDOM_TYPE FFTW_PATIENT

class ODEParams { // Class for ODE params
	private:
		
    public:
		int numT = -666;
		int numOmeg = -666;
		int doPlasmaCalc = 0;
		double chi_2, chi_3, mpi_sigmaK, mpi_k, ionE, sigmaBremsstrahlung, recombTime;
		double *omega, *kx, *rho, *y;
		complex<double> *k, *ee_p, *ee_m, *p_nl, *jhat, *j_e;
		fftw_plan p_ffft, j_ffft, ep_b, em_b, ep_f, em_f;

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
};




/*===========================================================================*/

class RootParams {
	private:
		int m_itnum = 1;
		int m_output = 0;
		int m_intCondition = 0; 
		int m_nRoot = -666;
		ODEParams *m_odeObj = nullptr;
    public:
		complex<double> *integral;
		RootParams() {} // Default constructor

      	RootParams(ODEParams *odeObj, int sizeRoot) { // Constructor Definition
			m_nRoot = sizeRoot;
			m_odeObj = odeObj;
			//m_odeObj = ODEParams(nT, omg);
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
		}

		/* ~RootParams() { // Destructor definition

			m_odeObj->~ODEParams();
		} */

		int getItNum() const { return m_itnum; }
		int getOutParam() const { return m_output; }
		int getIntCond() const { return m_intCondition; }
		int getSizeRoot() const { return m_nRoot; }
		ODEParams *getODEparams() const { return m_odeObj; }

		void setItNum(int aItNum) { m_itnum = aItNum; }
		void setOutParam(int aOutParam) { m_output = aOutParam; }
		void setIntCond(int aIntCond) { m_intCondition = aIntCond; }
		void setODEparams(ODEParams *aODEobj) { m_odeObj = aODEobj; }
};
