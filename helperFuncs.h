#include "BPPE.h"

void doLinearProblem(odeparam_type* p, complex<double> *sourceLeft, complex<double> *sourceRight, Structure& theStructure);
void generateTwoColorPulse(complex<double>* ee, fftw_plan e_f, complex<double>* source, pulseparam_type *pparams);
void boundary(double z, complex<double>*k_0, complex<double>*k_1, double *y);