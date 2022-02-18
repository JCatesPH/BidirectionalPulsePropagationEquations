#include "BPPE.h"

void doLinearProblem(ODEParams *odeObj, complex<double> *sourceLeft, complex<double> *sourceRight, Structure& theStructure);
void generateTwoColorPulse(complex<double>* ee, fftw_plan e_f, complex<double>* source, pulseparam_type *pparams);
void boundary(double z, complex<double>*k_0, complex<double>*k_1, double *y);
void generateGuess(gsl_vector *u, RootParams *rootObj, ODEParams *odeObj);
void updateGuess(double *ynew, complex<double> *sLeft, const gsl_vector *guessAm, RootParams *rootObj);
void dGdA(const gsl_vector *Am, void *rootparams, gsl_vector *df);