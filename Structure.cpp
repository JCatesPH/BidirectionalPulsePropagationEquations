//#include "stdafx.h"
#include "BPPE.h"
// Structure constructor

void boundary(double z, complex<double>*beta0, complex<double>*beta1, double *y) {

	complex<double> sp, sm;

    if (numDimensionsMinusOne == 1) {
        for (int i = 1; i < numOmX; i++)
        {
            complex<double> aPlus = y[i] + 1.0i*y[i + numOmX];
            complex<double> aMinus = y[i + 2 * numOmX] + 1.0i*y[i + 3 * numOmX];
            double kp2 = abs(kx[i] * conj(kx[i]));
            complex<double> denom = 2.0 * beta1[i] * (pow(beta1[i], 2) + kp2);

            sp = exp(1.0i*(beta0[i] - beta1[i])*z) * (beta0[i] + beta1[i]) * (beta0[i] * beta1[i] + kp2) / denom * aPlus 
                    + exp(-1.0i*(beta0[i] + beta1[i])*z) * (beta1[i] - beta0[i]) * (beta0[i] * beta1[i] - kp2) / denom * aMinus;

            sm = exp(1.0i*(beta0[i] + beta1[i])*z) * (beta1[i] - beta0[i]) * (beta0[i] * beta1[i] - kp2) / denom * aPlus 
                    + exp(-1.0i*(beta0[i] - beta1[i])*z) * (beta0[i] + beta1[i]) * (beta0[i] * beta1[i] + kp2) / denom * aMinus;
            
            y[i] = real(sp);
            y[i + numActiveOmega] = imag(sp);
            y[i + 2 * numActiveOmega] = real(sm);
            y[i + 3 * numActiveOmega] = imag(sm);
        }
    }   
    else {
        for (int i = 1; i < numActiveOmega; i++)
        {
            if (i >= freqLowerCutoff && i <= freqUpperCutoff) {
                complex<double> aPlus = y[i] + 1.0i*y[i + numActiveOmega];
                complex<double> aMinus = y[i + 2 * numActiveOmega] + 1.0i*y[i + 3 * numActiveOmega];

                sp = (exp(-1.0i*(beta0[i] - beta1[i])*z)*(beta0[i] + beta1[i]) / (2.0*beta1[i]) * aPlus + exp(1.0i*(beta0[i] + beta1[i])*z)*(beta1[i] - beta0[i]) / (2.0*beta1[i]) * aMinus);

                sm = (exp(-1.0i*(beta0[i] + beta1[i])*z)*(beta1[i] - beta0[i]) / (2.0*beta1[i]) * aPlus + exp(1.0i*(beta0[i] - beta1[i])*z)*(beta0[i] + beta1[i]) / (2.0*beta1[i]) * aMinus);
                
                y[i] = real(sp);
                y[i + numActiveOmega] = imag(sp);
                y[i + 2 * numActiveOmega] = real(sm);
                y[i + 3 * numActiveOmega] = imag(sm);
            }
        }
    }
	
	if (VERBOSE >= 5) { cout << "       Done Boundary() for z = " << z << endl;}
	return;
}
