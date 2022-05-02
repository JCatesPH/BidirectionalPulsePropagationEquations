#pragma once

#include "BPPE.h"

void writeInputEfield(std::complex<double>* ee_p);
void writeInputSpectrum(std::complex<double>* yp_init);
void write_out_eFieldAndSpectrumAtZlocation(int num, int j, double*y, double z, complex<double>*ee, complex<double>*k, fftw_plan e_b);
void write_multicolumnMonitor(int iterationNo, double theZpos, double *y, ODEParams *odeObj);
