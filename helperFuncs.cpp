#include "BPPE.h"

void doLinearProblem(ODEParams *odeObj, complex<double> *sourceLeft, complex<double> *sourceRight, Structure& theStructure) {
	for (int its = 0; its < 3; its++) {
		// Pass forward through boundaries 
		theStructure.doForwardPassThroughAllBoundaries(odeObj->y);

		// Set right source
		for (int i = 0; i < numActiveOmega; i++)
		{
			odeObj->y[i + 2*numActiveOmega] = real(sourceRight[i]);
			odeObj->y[i + 3*numActiveOmega] = imag(sourceRight[i]);
		}

		// Pass backward through boundaries
		theStructure.doBackwardPassThroughAllBoundaries(odeObj->y);

		// Reset left source to ensure consistency
		for (int i = 0; i < numActiveOmega; i++)
		{
			odeObj->y[i] = real(sourceLeft[i]);
			odeObj->y[i + numActiveOmega] = imag(sourceLeft[i]);
		}
	}
	
	// Output the reflected spectrum
	write_out_eFieldAndSpectrumAtZlocation(0, 0, odeObj->y, 0.0, odeObj->ee_m, theStructure.m_layers.front().getMaterial().getK(), odeObj->em_b);

	// Pass forward through boundaries 
	theStructure.doForwardPassThroughAllBoundaries(odeObj->y);

	// Output the transmitted spectrum 
	write_out_eFieldAndSpectrumAtZlocation(0, 1, odeObj->y, theStructure.getThickness(), odeObj->ee_p, theStructure.m_layers.back().getMaterial().getK(), odeObj->ep_b);

	// Set right source
	/* for (int i = 0; i < numActiveOmega; i++)
	{
		y[i + 2*numActiveOmega] = real(ym_init[i]);
		y[i + 3*numActiveOmega] = imag(ym_init[i]);
	} */

	// Pass backward through boundaries
	theStructure.doBackwardPassThroughAllBoundaries(odeObj->y);

	// Reset left source to ensure consistency
	for (int i = 0; i < numActiveOmega; i++)
	{
		odeObj->y[i] = real(sourceLeft[i]);
		odeObj->y[i + numActiveOmega] = imag(sourceLeft[i]);
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
	if(fetestexcept(FE_OVERFLOW | FE_INVALID | FE_DIVBYZERO)) raise(SIGFPE);
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


void generateGuess(gsl_vector *u, RootParams *rootObj, ODEParams *odeObj) {
	int sizeRoot = rootObj->getSizeRoot();
	int itSecant = 1;
	//ODEParams *odeObj = rootObj->getODEparams();
	gsl_vector *tmp = gsl_vector_alloc(sizeRoot);
	complex<double> map1, map2, initGuess;

	// Create pseudo-random number generator for initial guess
	/* random_device rd;
	mt19937 gen(rd()); */
	default_random_engine gen;
	uniform_real_distribution<double> dis(-NOISE_MAGNITUDE, NOISE_MAGNITUDE);
	
	/* complex<double>* integral1 = (complex<double>*)malloc(sizeof(complex<double>) * numActiveOmega);
	complex<double>* integral2 = (complex<double>*)malloc(sizeof(complex<double>) * numActiveOmega);
	complex<double>* Am_guess1 = (complex<double>*)malloc(sizeof(complex<double>) * numActiveOmega);
	complex<double>* Am_guess2 = (complex<double>*)malloc(sizeof(complex<double>) * numActiveOmega);

	for (int i = 0; i < numActiveOmega; i++)
	{
		integral1[i] = 0.0;
		integral2[i] = 0.0;
	}

	// Use integral condition to inform guess
    rootObj->setOutParam(1); // Turn output on (1) or off (0)
    rootObj->setIntCond(1); // Set whether to calculate integral
	rootObj->integral = integral1;
	mapG(u, rootObj, tmp);

	for (int k = 0; k < numActiveOmega; k++) {
		Am_guess1[k] = odeObj->y[k + 2 * numActiveOmega] + 1.0i * odeObj->y[k + 3 * numActiveOmega];
	}


	// Get second guess for secant method
    rootObj->setItNum(2);
    rootObj->setIntCond(2);
	rootObj->integral = integral2;
	mapG(u, rootObj, tmp);

	for (int k = 0; k < numActiveOmega; k++) {
		Am_guess2[k] = odeObj->y[k + 2 * numActiveOmega] + 1.0i * odeObj->y[k + 3 * numActiveOmega];
	}

	// Do Secant update
	for (int k = 0; k < sizeRoot/2; k++){
		map1 = Am_guess1[k + freqLowerCutoff] + integral1[k + freqLowerCutoff];
		map2 = Am_guess2[k + freqLowerCutoff] + integral2[k + freqLowerCutoff];
		if (abs(map2 - map1) < DBL_MIN) {
			gsl_vector_set(u, k, real(Am_guess2[k + freqLowerCutoff]));
			gsl_vector_set(u, k + sizeRoot/2, imag(Am_guess2[k + freqLowerCutoff]));
		}
		else {
			initGuess = (Am_guess1[k + freqLowerCutoff] * map2 - Am_guess2[k + freqLowerCutoff] * map1) 
			/ (map2 - map1);
			gsl_vector_set(u, k, real(initGuess));
			gsl_vector_set(u, k + sizeRoot/2, imag(initGuess));
		}
	}
	//cout << "Finished iteration 3" << endl;
	//if(fetestexcept(FE_OVERFLOW | FE_INVALID | FE_DIVBYZERO)) raise(SIGFPE); 

	for (int it = 3; it < itSecant; it++) {
		// Move latest update to previous guess
		for (int k = 0; k < numActiveOmega; k++) {
			Am_guess1[k] = Am_guess2[k]; 
			integral1[k] = integral2[k];
			integral2[k] = 0.0;
		}

		// Do another iteration
		rootObj->setItNum(it);
		rootObj->setIntCond(it); // Set whether to calculate integral
		rootObj->integral = integral2;
		mapG(u, rootObj, tmp);

		for (int k = 0; k < numActiveOmega; k++) {
			Am_guess2[k] = odeObj->y[k + 2 * numActiveOmega] + 1.0i * odeObj->y[k + 3 * numActiveOmega];
		}

		// Do Secant update
		for (int k = 0; k < sizeRoot/2; k++){
			map1 = Am_guess1[k + freqLowerCutoff] + integral1[k + freqLowerCutoff];
			map2 = Am_guess2[k + freqLowerCutoff] + integral2[k + freqLowerCutoff];
			if (abs(map2 - map1) < DBL_MIN) {
				gsl_vector_set(u, k, real(Am_guess2[k + freqLowerCutoff]));
				gsl_vector_set(u, k + sizeRoot/2, imag(Am_guess2[k + freqLowerCutoff]));
			}
			else {
				initGuess = (Am_guess1[k + freqLowerCutoff] * map2 - Am_guess2[k + freqLowerCutoff] * map1) 
				/ (map2 - map1);
				gsl_vector_set(u, k, real(initGuess));
				gsl_vector_set(u, k + sizeRoot/2, imag(initGuess));
			}
		}
	} 

	// Free some memory
	free(Am_guess1);
	free(Am_guess2);
	free(integral1);
	free(integral2);
	gsl_vector_free(tmp); */
	
	
	/* cout << "Finished iteration 4" << endl;
	if(fetestexcept(FE_OVERFLOW | FE_INVALID | FE_DIVBYZERO)) raise(SIGFPE); */
	for (int k = 0; k < sizeRoot/2; k++){
		//gsl_vector_set(u, k, (1.0 + dis(gen)) * odeObj->y[k + 2*numActiveOmega + freqLowerCutoff]);
		//gsl_vector_set(u, k + sizeRoot/2, (1.0 + dis(gen)) * odeObj->y[k + 3*numActiveOmega + freqLowerCutoff]);
		//gsl_vector_set(u, k, (1.0 + dis(gen)) * gsl_vector_get(u, k));
		//gsl_vector_set(u, k, (1.0 + dis(gen)) * gsl_vector_get(u, k + sizeRoot/2));
		gsl_vector_set(u, k, dis(gen));
		gsl_vector_set(u, k + sizeRoot/2, dis(gen));
	}


	rootObj->setItNum(itSecant);
	rootObj->setOutParam(0); // Turn output on (1) or off (0)
    rootObj->setIntCond(0); // Set whether to calculate integral

}



void updateGuess(double *ynew, complex<double> *sLeft, const gsl_vector *guessAm, RootParams *rootObj) {
	for (int k = 0; k < numActiveOmega; k++){
		ynew[k] = real(sLeft[k]);
		ynew[k + numActiveOmega] = imag(sLeft[k]);
	}
	for (int k = 0; k < rootObj->getSizeRoot()/2; k++){
		ynew[k + 2*numActiveOmega + freqLowerCutoff] = gsl_vector_get(guessAm, k);
		ynew[k + 3*numActiveOmega + freqLowerCutoff] = gsl_vector_get(guessAm, k + rootObj->getSizeRoot()/2);
	}

}

void dGdA(const gsl_vector *Am_in, void *rootparams, gsl_vector *dG) {
	RootParams *rootObj = reinterpret_cast<RootParams*>(rootparams);
	const int nRoot = rootObj->getSizeRoot();
	//double epsrel = 1e-4;
	double epsrel = GSL_SQRT_DBL_EPSILON; // 1.49e-8

	const double G0 = rootObj->getGnorm();
	gsl_vector *Am_1 = gsl_vector_alloc(nRoot);
	gsl_vector_memcpy(Am_1, Am_in);

	for (int j = 0; j < nRoot; j++){
		double Aj = gsl_vector_get(Am_1, j);
		double deltaA = epsrel * abs(Aj);
		//double deltaA = epsrel;

		if (deltaA == 0) {
			deltaA = epsrel;
		}

		gsl_vector_set(Am_1, j, Aj + deltaA);

		double G1 = mapG(Am_1, rootparams);

		//gsl_vector_set(dG, j, (G1-G0)/(deltaA + 1e-5));
		gsl_vector_set(dG, j, (G1-G0)/(deltaA));
	}
	
}

void GdG(const gsl_vector *Am_in, void *rootparams, double *G, gsl_vector *dG) {
	*G = mapG(Am_in, rootparams);
	dGdA(Am_in, rootparams, dG);
}
