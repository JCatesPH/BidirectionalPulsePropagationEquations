#pragma once

#include "BPPE.h"

void doLinearProblem(ODEParams *odeObj, complex<double> *sourceLeft, complex<double> *sourceRight, Structure& theStructure);
void generateTwoColorPulse(complex<double>* ee, fftw_plan e_f, complex<double>* source, pulseparam_type *pparams);
void generateGuess(gsl_vector *u, RootParams *rootObj, ODEParams *odeObj);
void updateGuess(double *ynew, complex<double> *sLeft, const gsl_vector *guessAm, RootParams *rootObj);
