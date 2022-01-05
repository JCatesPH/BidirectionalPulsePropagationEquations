#include "BPPE.h"

void doLinearProblem(odeparam_type* p, complex<double> *sourceLeft, complex<double> *sourceRight, Structure& theStructure) {
	for (int its = 0; its < 3; its++) {
		// Pass forward through boundaries 
		theStructure.doForwardPassThroughAllBoundaries(p->y);

		// Set right source
		for (int i = 0; i < numActiveOmega; i++)
		{
			p->y[i + 2*numActiveOmega] = real(sourceRight[i]);
			p->y[i + 3*numActiveOmega] = imag(sourceRight[i]);
		}

		// Pass backward through boundaries
		theStructure.doBackwardPassThroughAllBoundaries(p->y);

		// Reset left source to ensure consistency
		for (int i = 0; i < numActiveOmega; i++)
		{
			p->y[i] = real(sourceLeft[i]);
			p->y[i + numActiveOmega] = imag(sourceLeft[i]);
		}
	}
	
	// Output the reflected spectrum
	write_out_eFieldAndSpectrumAtZlocation(0, 0, p->y, 0.0, p->ee_m, theStructure.m_layers.front().getMaterial().getK(), p->em_b);

	// Pass forward through boundaries 
	theStructure.doForwardPassThroughAllBoundaries(p->y);

	// Output the transmitted spectrum 
	write_out_eFieldAndSpectrumAtZlocation(0, 1, p->y, theStructure.getThickness(), p->ee_p, theStructure.m_layers.back().getMaterial().getK(), p->ep_b);

	// Set right source
	/* for (int i = 0; i < numActiveOmega; i++)
	{
		y[i + 2*numActiveOmega] = real(ym_init[i]);
		y[i + 3*numActiveOmega] = imag(ym_init[i]);
	} */

	// Pass backward through boundaries
	theStructure.doBackwardPassThroughAllBoundaries(p->y);

	// Reset left source to ensure consistency
	for (int i = 0; i < numActiveOmega; i++)
	{
		p->y[i] = real(sourceLeft[i]);
		p->y[i + numActiveOmega] = imag(sourceLeft[i]);
	}

}

void generateTwoColorPulse(complex<double>* ee, fftw_plan e_f, complex<double>* source, pulseparam_type *pparams) {
    double ht = (1.0 * domain_t) / double(num_t);
	int o;

    for (int i = 0; i < num_t; i++) {
        double s = (-domain_t/2.0) + ht * ((double)i + 1);
        ee[i] = pparams->A0 * (sqrt(1.0 - pparams->relativeIntensity) * exp(-2.0 * log(2.0) * pow(s, 2) / (pow(pparams->pulseDuration, 2))) * cos(pparams->omega0 * s)
                        + sqrt(pparams->relativeIntensity) *       exp(-8.0 * log(2.0) * pow(s, 2) / (pow(pparams->pulseDuration, 2))) * cos(2.0 * pparams->omega0 * s + pparams->relativePhase));
		//ee[i] = 0.0;
	}
    writeInputEfield(ee);

    fftw_execute(e_f);
    normalizeFFT(ee);

	for (int i = 0; i < numActiveOmega; i++) {
		source[i] = ee[i];
	}

}


void boundary(double z, complex<double>*k_0, complex<double>*k_1, double *y) {

	complex<double> sp, sm;
	for (int i = 1; i < numActiveOmega; i++)
	{
		if (i >= freqLowerCutoff && i <= freqUpperCutoff) {
			complex<double> aPlus = y[i] + 1.0i*y[i + numActiveOmega];
			complex<double> aMinus = y[i + 2 * numActiveOmega] + 1.0i*y[i + 3 * numActiveOmega];

			sp = (exp(1.0i*(k_0[i] - k_1[i])*z)*(k_0[i] + k_1[i]) / (2.0*k_1[i]) * aPlus + exp(-1.0i*(k_0[i] + k_1[i])*z)*(k_1[i] - k_0[i]) / (2.0*k_1[i]) * aMinus);

			sm = (exp(1.0i*(k_0[i] + k_1[i])*z)*(k_1[i] - k_0[i]) / (2.0*k_1[i]) * aPlus + exp(-1.0i*(k_0[i] - k_1[i])*z)*(k_0[i] + k_1[i]) / (2.0*k_1[i]) * aMinus);
			
			y[i] = real(sp);
			y[i + numActiveOmega] = imag(sp);
			y[i + 2 * numActiveOmega] = real(sm);
			y[i + 3 * numActiveOmega] = imag(sm);
		}
	}
	if (VERBOSE >= 5) { cout << "       Done Boundary() for z = " << z << endl;}
	return;
}