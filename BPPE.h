#pragma once

//#include "stdafx.h"
#define _USE_MATH_DEFINES

#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <iomanip> 
#include <complex> 
#include <cmath>
#include <algorithm>
#include <omp.h>
#include <fftw3.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_vector.h>
#include "gsl/gsl_multiroots.h"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include "gsl/gsl_math.h"
#include <string>
#include "physicalConstants.h"
#include "Structure.h"
#include "Utilities.h"
#include "gslStructs.h"

using namespace std;
#define STRING_BUFFER_SIZE 256

//#define FFTW_WISDOM_TYPE FFTW_ESTIMATE
#define FFTW_WISDOM_TYPE FFTW_PATIENT
#define DO_DRUDE_MODEL
//#define DO_CONSTPLASMA
//#define DO_ARGON_PLASMA

#define NOISE_MAGNITUDE 1.0e-6

// CODE parameters
#define USE_CPP_BOUNDARY
#define USE_CPP_NONLINEAR
//#define FIXEDSTEP

#undef  WRITE_OUT_REFLECTANCE

#ifdef WRITE_OUT_REFLECTANCE
complex<double>* eFieldPlusBACKUPCOLM;
#endif


// Simulation parameters
extern char SIM_DATA_OUTPUT[30];
extern int VERBOSE;
const int num_Threads = 4; // numnber of OpenMP threads
const int num_iterations = 5; //number of BPPE iterations
const int numDimensionsMinusOne = 0; //(1+1) dimension (0) or (2+1) dimension (1)
const int normType = -1;
//const double zStepMaterial1 = 0.01 * microns; // Step for integration
extern double zStepMaterial1;

// GSL ODE API parameters
const double ode_epsabs = 1e-9;
const double ode_epsrel = 1e-6;
const int ode_nmax = 1e6;

// GSL Quasi-Newton API parameters
const double root_epsabs = 1e-9;
const double root_epsrel = 1e-6;


// plasma parameters
//const int plasmaOnOff = 0; //plasma off (0) Using Andrew (1) Using UPPE MPI (2)
const double num_atoms = 2.0e25;  //number of atoms in gas [1/m^3]
const double rho_0 = 9.0e24; // initial electron density
const double j_e0 = 0.0;
const double omegaPlasmaDamping = 2.0 * M_PI * 5.3e14; //2.0 * M_PI*5.3e12;  //plasma damping
const double tauCollision = 190.0e-15; // <- This number used in Berge paper. //26.9984566e-15; //1.88679e15; //mean collision time
const double mpi_sigmaK = 3.4e-128;
const double mpi_k = 8.0;
const double FUDGE_FACTOR = 1.0;
const double Znaught = FUDGE_FACTOR * (1.0 / (epsilon_0 * cLight));
//const double nu_a = 4.13e16; // [Hz]
//const double E_a = 5.14e11; // [V/m]
//const double U_H = 13.6; // Ionization potential of hydrogen [eV]


// domain parameters
extern int num_t;
extern double domain_t;
extern int numActiveOmega;
const int num_x = int(pow(2, 6)); //number of x points
const double domain_x = 125.0e-6; //x domain
extern int freqLowerCutoff; 
extern int freqUpperCutoff;
const double shift = 1.0e-9;
//const int numActiveOmega = num_t / 2 + 1; //num_x * num_t / 2 + 2;
//const int numActiveOmega2 = numActiveOmega - (num_t / 2 + 1);
//const int l_0 = (num_t / 2 + 1)*(num_x / 2 + 1);

extern double LHSsourceLayerThickness; //distance from laser source to slab
extern double RHSbufferLayerThickness; //distance from slab to receiver
extern double sampleLayerThickness;

// material parameters

// Predefine a few common material parameters
const double n0_Vacuum = 1.0; //central index in material 0
const double n0_Argon = 1.00026436; // For 1.0 um, from https://refractiveindex.info/?shelf=main&book=Ar&page=Peck-15C
const double n2_Argon = 5.0e-23; //nonlinear index in material 1; ~1e-19 cm^2/W from doi: 10.1007/s00340-013-5354-0
const double chi3_Argon = (4 / 3) * epsilon_0 * cLight * pow(n0_Argon, 2) * n2_Argon;
const double chi2_Argon = 0.0;
const double n0_Material1 = 2.5; //1.5; //central index in material 1
const double n0_Material2 = 1.5; //1.5; //central index in material 2
const double n0_Material3 = 1.0; //central index in material 3
const double n2_Material1 = 0.0; //nonlinear index in material 1
const double n2_Material2 = 0.0; //nonlinear index in material 2
const double chi_2 = 0.0; //30.0e-12; //nonlinear chi_2
const double chi3_Material1 = (4 / 3) * epsilon_0 * cLight * pow(n0_Material1, 2) * n2_Material1;
const double chi3_Material2 = (4 / 3) * epsilon_0 * cLight * pow(n0_Material2, 2) * n2_Material2;
const double chi2_Material1 = 0.0;
const double chi2_Material2 = 0.0; //-1.0 * chi_2;

//const double A_0 = sqrt(2.0 * I_0 / (epsilon_0*n0_Vacuum*cLight)); //pulse peak amplitude
//const double Keldysh = omega_0 * sqrt(2.0 * u_Argon*charge_e) / A_0; //keldysh parameter
const double omegaPlasma = sqrt(num_atoms*pow(charge_e,2)/(mass_e*epsilon_0)); //fully-ionized plasma frequency

//sellemeier parameters
const double Sellmeir_chi_1_1 = 2.4272;
const double Sellmeir_chi_1_2 = 1.4617;
const double Sellmeir_chi_1_3 = 9.6536;
const double Sellmeir_omega_1 = 1.5494e16;
const double Sellmeir_omega_2 = 7.9514e15;
const double Sellmeir_omega_3 = 9.7766e13;


typedef struct {
	double A0;
	double relativeIntensity;
	double relativePhase;
	double pulseDuration;
	double omega0;
}pulseparam_type;


//void doNonlinearPartofBPPE();
void iterateBPPE();
void writeInputEfield(std::complex<double>* ee_p);
void writeInputSpectrum(std::complex<double>* yp_init);
//void generateTwoColorPulse(complex<double>* ee, fftw_plan e_f, complex<double>* source, pulseparam_type *pparams);
//void initalizeYarray(double* y, std::complex<double>* yp_init, std::complex<double>* ym0_init);
//void write2DtoFile(std::complex<double>* ee_p);
//void fillYfromYpAndYm(double* y, std::complex<double>* yp_init, std::complex<double>* ym0_init);
//void initalizeArrays(std::complex<double>* ym1_init, std::complex<double>* ym0_init, std::complex<double>* integral);
void boundary(double z, complex<double>*k_0, complex<double>*k_1, double *y);
//void update_guess(complex<double>*yp_init, complex<double>*f0, complex<double>*ym0_init, double*y, complex<double>*integral);
void writeSimParameters();
void generateLayers(MaterialDB &myMaterialsDB, Structure &theStructure);
void setupPointMonitorLocations(MaterialDB& myMaterialsDB, Structure& theStructure);
//void set_guess(complex<double>* ee_p, complex<double>* yp_init, complex<double>* ym0_init, complex<double>* ym1_init, complex<double>* ym1_temp, double* y, fftw_plan ep_f, complex<double>* ee_m, fftw_plan em_b, fftw_plan ep_b, complex<double>* integral);
void write_out_eFieldAndSpectrumAtZlocation(int num, int j, double*y, double z, complex<double>*ee, complex<double>*k, fftw_plan e_b);
//void am_to_zero(double*y);
//int new_initial_data(complex<double>*ym0_init, complex<double>*ym1_init, complex<double>*ym1_temp, complex<double>*yp_init, double*y, complex<double>*integral);
int func(double z, const double y[], double f[], void *params);
void integrate(double z, double zStep, odeparam_type *params, double*y, complex<double>*integral);
//void write_multicolumnMonitor(int iterationNo, double theZpos, complex<double>* eep, complex<double>* eem, double* ne, complex<double>* j_e);
void write_multicolumnMonitor(int iterationNo, double theZpos, double *y, odeparam_type *p);
void DELME_ArgonDispersion(double* omg);
void DELME_AndrewPreformed(double* omg, Material* mat);
