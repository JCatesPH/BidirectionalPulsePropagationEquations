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
#include "omp.h"
#include <fftw3.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_math.h>
#include <string>
#include"physicalConstants.h"
#include"Structure.h"

using namespace std;
#define STRING_BUFFER_SIZE 256
#define INITIAL_GUESS_SEED_VALUE 1.0e5

// CODE parameters
int VERBOSE = 6;
#define USE_CPP_BOUNDARY
#define USE_CPP_NONLINEAR

#undef  WRITE_OUT_REFLECTANCE

#ifdef WRITE_OUT_REFLECTANCE
complex<double>* eFieldPlusBACKUPCOLM;
#endif


// Simulation parameters
#ifdef _WIN32
const char *SIM_DATA_OUTPUT=".\\DATA\\";
#endif

#ifdef __unix__
const char *SIM_DATA_OUTPUT="./DATA/";
#endif

//pulse parameters
const double I_0 = 1.5e17;  //initial peak intensity
const double twoColorSH_amplitude = 0.1;  //0.1; //two-color pulse: 2nd harmonic with half duration of fundamental
const double twoColorSH_phase = M_PI_2; //phase shift of 2nd harmonic
const double tau = 50.0e-15; //pulse duration fwhm
const double lambda_0 = 0.8e-6; //central wavelength
const double omega_0 = 2 * M_PI*cLight / lambda_0; //central angular frequency
const double waist_x = 50.0e-6; //pulse-waist fwhm
const int num_t = int(pow(2, 14)); //number of time points
const int num_x = int(pow(2, 6)); //number of x points
const double domain_t = 1500e-15; //time domain
const double domain_x = 125.0e-6; //x domain
const double num_atoms = 4.0e25;  //number of atoms in gas
const double rho_0 = 1000; // = num_atoms; // initial electron density
const double j_e0 = 0.0;
const int freqUpperCutoff = num_t / 2; // upper frequency cut-off
const int freqLowerCutoff = 1; //lower frequency cut-off
const int num_iterations = 5; //number of BPPE iterations
const int num_Threads = 4; // numnber of OpenMP threads
const double shift = 1.0e-9;

//material parameters
const int numDimensionsMinusOne = 0; //(1+1) dimension (0) or (2+1) dimension (1)
const int plasmaOnOff = 2; //plasma off (0) Using Andrew (1) Using UPPE MPI (2)
const int l_0 = (num_t / 2 + 1)*(num_x / 2 + 1);
const int numActiveOmega = num_t / 2 + 1; //num_x * num_t / 2 + 2;
const int numActiveOmega2 = numActiveOmega - (num_t / 2 + 1);
const double omegaPlasmaDamping = 2.0 * M_PI*5.3e14; //2.0 * M_PI*5.3e12;  //plasma damping
const double tauCollision = 26.9984566e-15; //1.88679e15; //mean collision time
const double lengthSample = 5.0e-6;  //length of slab
const double distanceSourceToSample = 10.0e-6; //distance from laser source to slab
const double distanceSampleToReceiver = 10.0e-6; //distance from slab to receiver
const double zRightHandSideOfSample = distanceSourceToSample+lengthSample;
const double Z_4 = zRightHandSideOfSample; 
const double sampleLayerThickness = 0.25e-6; //half period
const int numLayersInSample = int(lengthSample / sampleLayerThickness); //number of half-periods in slab
const double zStepMaterial1 = 0.01e-6; //propagation step in material 1
const double zStepMaterial2 = 0.01e-6; //propagation step in material 2
const int numZstepsMaterial1 = int(sampleLayerThickness / zStepMaterial1); //number of steps in material 1
const int numZstepsMaterial2 = int(sampleLayerThickness / zStepMaterial2); //number of steps in material 2

const double n0_Material0 = 1.0; //central index in material 0
const double n0_Material1 = 1.750; //1.5; //central index in material 1
const double n0_Material2 = 1.25; //1.5; //central index in material 2
const double n0_Material3 = 1.0; //central index in material 3
const double n2_Material1 = 2.0e-20; //nonlinear index in material 1
const double n2_Material2 = 4.0e-20; //nonlinear index in material 2
const double chi_2 = 0.0; //30.0e-12; //nonlinear chi_2
const double A_0 = sqrt(2.0 * I_0 / (epsilon_0*n0_Material0*cLight)); //pulse peak amplitude
const double Keldysh = omega_0 * sqrt(2.0 * u_Argon*charge_e) / A_0; //keldysh parameter
const double omegaPlasma = sqrt(num_atoms*pow(charge_e,2)/(mass_e*epsilon_0)); //fully-ionized plasma frequency

//plasma parameters
const double mpi_sigmaK = 3.4e-128;
const double mpi_k = 8.0;
const double FUDGE_FACTOR = 1.0;
const double Znaught = FUDGE_FACTOR * (1.0 / (epsilon_0 * cLight));

//sellemeier parameters
const double Sellmeir_chi_1_1 = 2.4272;
const double Sellmeir_chi_1_2 = 1.4617;
const double Sellmeir_chi_1_3 = 9.6536;
const double Sellmeir_omega_1 = 1.5494e16;
const double Sellmeir_omega_2 = 7.9514e15;
const double Sellmeir_omega_3 = 9.7766e13;

typedef struct {

	double chi_2;
	double chi_3;
	double *omega, *kx, *rho;
	complex<double> *k, *ee_p, *ee_m, *nl_k, *nl_p, *j_e;
	fftw_plan nk_f, ep_b, em_b, np_f;

}param_type;

complex<double> index_0(double omg, double kx);
void doNonlinearPartofBPPE();
complex<double> index_1(double omg, double kx, errno_t err, FILE*fp);
complex<double> index_2(double omg, double kx);
complex<double> index_3(double omg, double kx);
param_type* fill_params(double chi_2, double chi_3, double*omg, double*kx, double*ne, complex<double>*j_e, complex<double>*k, complex<double>*ee_p, complex<double>*ee_m, complex<double>*nl_k, complex<double>*nl_p, fftw_plan nk_f, fftw_plan ep_b, fftw_plan em_b, fftw_plan np_f);
complex<double> exact_linear(complex<double> *k_0, complex<double> *k_1, complex<double> *k_2, complex<double> *k_3, complex<double>*k_4, int i);
void fill_omg_k(double*omg, double*kx, complex<double>*k_0, complex<double>*k_1, complex<double>*k_2, complex<double>*k_3, errno_t err, FILE*fp);
void writeInputEfield(std::complex<double>* ee_p);
void writeInputSpectrum(std::complex<double>* yp_init);
void initalizeYarray(double* y, std::complex<double>* yp_init, std::complex<double>* ym0_init);
void write2DtoFile(std::complex<double>* ee_p);
void fillYfromYpAndYm(double* y, std::complex<double>* yp_init, std::complex<double>* ym0_init);
void set_guess(complex<double>*ee_p, complex<double>*yp_init, complex<double>*ym0_init, complex<double>*ym1_init, complex<double>*ym1_temp, complex<double>*f0, complex<double>*f1, double*y, fftw_plan ep_f, complex<double>*ee_m, fftw_plan em_b, fftw_plan ep_b, complex<double>*k_0, complex<double>*k_1, complex<double>*k_2, complex<double>*k_3, complex<double>*integral);
void initalizeArrays(std::complex<double>* ym1_init, std::complex<double>* ym0_init, std::complex<double>* integral);
void boundary(double z, complex<double>*k_0, complex<double>*k_1, double *y);
void update_guess(complex<double>*yp_init, complex<double>*f0, complex<double>*ym0_init, double*y, complex<double>*integral);
void writeSimParameters();
void generateApp1MaterialsAndStructure(MaterialDB &myMaterialsDB, Structure &theStructure);
void generateDefectMaterialsAndStructure(MaterialDB& myMaterialsDB, Structure& theStructure);
void generatePlasmaTestMaterialsAndStructure(MaterialDB& myMaterialsDB, Structure& theStructure);

void write_out_eFieldAndSpectrumAtZlocation(int num, int j, double*y, double z, complex<double>*ee, complex<double>*k, fftw_plan e_b);
void am_to_zero(double*y);
int new_initial_data(complex<double>*ym0_init, complex<double>*ym1_init, complex<double>*ym1_temp, complex<double>*yp_init, complex<double>*f0, complex<double>*f1, double*y, complex<double>*integral);
int func(double z, const double y[], double f[], void *params);
void integrate(double z, double chi_2, double chi_3, complex<double>*k, double*omg, double*ne, complex<double>*j_e, double*y, complex<double>*integral, complex<double>*ee_p, complex<double>*ee_m, complex<double>*nl_k, complex<double>*nl_p, fftw_plan ep_b, fftw_plan em_b, fftw_plan nk_f, fftw_plan np_f);
void write_multicolumnMonitor(int iterationNo, double theZpos, complex<double>* eep, complex<double>* eem, double* ne, complex<double>* j_e);
void write_out_ne(int j, double zPos, double* ne);
void write_out_je(int j, complex<double>*j_e);
void write_out_ee_p(int j, complex<double>* eep);
void write_out_ee_m(int j, complex<double>* eem);
