#include "BPPE.h"

void doLinearProblem(ODEParams *odeObj, complex<double> *sourceLeft, complex<double> *sourceRight, Structure& theStructure) {
	complex<double>* Am1 = (complex<double>*)malloc(sizeof(complex<double>) * numActiveOmega);
	complex<double>* Am2 = (complex<double>*)malloc(sizeof(complex<double>) * numActiveOmega);
	complex<double>* GAm1 = (complex<double>*)malloc(sizeof(complex<double>) * numActiveOmega);
	complex<double>* GAm2 = (complex<double>*)malloc(sizeof(complex<double>) * numActiveOmega);

	write_out_eFieldAndSpectrumAtZlocation(0, 0, odeObj->y, 0.0, odeObj->ee_m, theStructure.m_layers.front().getMaterial().getK(), odeObj->em_b);
	// Pass forward through boundaries 
	theStructure.doForwardPassThroughAllBoundaries(odeObj->y);
	write_out_eFieldAndSpectrumAtZlocation(0, 1, odeObj->y, theStructure.getThickness(), odeObj->ee_p, theStructure.m_layers.back().getMaterial().getK(), odeObj->ep_b);

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

	// -------------------------------
	// Set first guess to LHS Am
	for (int k = 0; k < numActiveOmega; k++) {
		Am1[k] = odeObj->y[k + 2 * numActiveOmega] + 1.0i * odeObj->y[k + 3 * numActiveOmega];
	}

	// Pass forward through boundaries 
	theStructure.doForwardPassThroughAllBoundaries(odeObj->y);

	// Set map of LHS Am to RHS
	for (int j = 0; j < numActiveOmega; j++) {
		GAm1[j] = (odeObj->y[j + 2 * numActiveOmega] + 1.0i * odeObj->y[j + 3 * numActiveOmega]) * exp(-1.0i * odeObj->k[j] * theStructure.getThickness());
	}

	theStructure.doBackwardPassThroughAllBoundaries(odeObj->y);

	// -------------------------------
	// Set second guess with r.v.
	default_random_engine gen;
	uniform_real_distribution<double> dis(0.0, 1e2);

	// Reset left source to ensure consistency
	for (int i = 0; i < numActiveOmega; i++)
	{
		Am2[i] = dis(gen) + 1.0i * dis(gen);

		odeObj->y[i] = real(sourceLeft[i]);
		odeObj->y[i + numActiveOmega] = imag(sourceLeft[i]);
		odeObj->y[i + 2*numActiveOmega] = real(Am2[i]);
		odeObj->y[i + 3*numActiveOmega] = imag(Am2[i]);
	}

	// Set map of LHS Am to RHS
	for (int j = 0; j < numActiveOmega; j++) {
		GAm2[j] = (odeObj->y[j + 2 * numActiveOmega] + 1.0i * odeObj->y[j + 3 * numActiveOmega]) * exp(-1.0i * odeObj->k[j] * theStructure.getThickness());
	}

	theStructure.doBackwardPassThroughAllBoundaries(odeObj->y);

	// -------------------------------
	for (int its = 0; its < 3; its++) {
		for (int k = 0; k < numActiveOmega; k++){
			odeObj->y[k] = real(sourceLeft[k]);
			odeObj->y[k + numActiveOmega] = imag(sourceLeft[k]);

			complex<double> secant;
			if(abs(GAm1[k] - GAm2[k]) < DBL_MIN) {
				secant = Am2[k];
			}

			else{
				secant = (Am1[k] * GAm2[k] - Am2[k] * GAm1[k]) / (GAm2[k] - GAm1[k]);
			}

			odeObj->y[k + 2*numActiveOmega] = real(secant);
			odeObj->y[k + 3*numActiveOmega] = imag(secant);

			Am1[k] = Am2[k];
			GAm1[k] = GAm2[k];
			Am2[k] = secant;
		}
		// Pass forward through boundaries 
		theStructure.doForwardPassThroughAllBoundaries(odeObj->y);

		// Set map of LHS Am to RHS
		for (int j = 0; j < numActiveOmega; j++) {
			GAm2[j] = (odeObj->y[j + 2 * numActiveOmega] + 1.0i * odeObj->y[j + 3 * numActiveOmega]) * exp(-1.0i * odeObj->k[j] * theStructure.getThickness());
		}

		// Set right source
		for (int i = 0; i < numActiveOmega; i++)
		{
			odeObj->y[i + 2*numActiveOmega] = real(sourceRight[i]);
			odeObj->y[i + 3*numActiveOmega] = imag(sourceRight[i]);
		}

		// Pass backward through boundaries
		theStructure.doBackwardPassThroughAllBoundaries(odeObj->y);
	}
	
	// Output the reflected spectrum
	//write_out_eFieldAndSpectrumAtZlocation(0, 0, odeObj->y, 0.0, odeObj->ee_m, theStructure.m_layers.front().getMaterial().getK(), odeObj->em_b);
	
	// Reset left source to ensure consistency
	for (int i = 0; i < numActiveOmega; i++)
	{
		odeObj->y[i] = real(sourceLeft[i]);
		odeObj->y[i + numActiveOmega] = imag(sourceLeft[i]);
	}

	// Pass forward through boundaries 
	theStructure.doForwardPassThroughAllBoundaries(odeObj->y);

	// Output the transmitted spectrum 
	//write_out_eFieldAndSpectrumAtZlocation(0, 1, odeObj->y, theStructure.getThickness(), odeObj->ee_p, theStructure.m_layers.back().getMaterial().getK(), odeObj->ep_b);

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


void generateTwoColorPulse(complex<double>* ee, fftw_plan e_f, complex<double>* source, pulseparam_type *pparams) {
    double ht = domain_t / double(num_t);
	int o;
	if (numDimensionsMinusOne == 1) {
		double hx = domain_x / double(num_x);
		for (int j = 0; j < num_x; j++) {
			for (int i = 0; i < num_t; i++) {
				double r = (-domain_x/2.0) + hx * ((double)j + 1);
				double s = (-domain_t/2.0) + ht * ((double)i + 1);
				ee[i + num_t * j] = pparams->A0 * (sqrt(1.0 - pparams->relativeIntensity) * exp(-2.0 * log(2.0) * pow(s, 2) / pow(pparams->pulseDuration, 2)) * exp(-2.0 * log(2.0) * pow(r, 2) / pow(pparams->pulseWaist, 2)) * cos(pparams->omega0 * s)
								+ sqrt(pparams->relativeIntensity) * exp(-8.0 * log(2.0) * pow(s, 2) / (pow(pparams->pulseDuration, 2))) * exp(-2.0 * log(2.0) * pow(r, 2) / pow(pparams->pulseWaist, 2)) * cos(2.0 * pparams->omega0 * s + pparams->relativePhase));
				}
		}
		writeInputEfield(ee);

		fftw_execute(e_f);
		normalizeFFT(ee, fftnorm);
		
		for (int j = 0; j <= num_x / 2; j++) {
			if (j == 0) {
				for (int i = 0; i < numActiveOmega; i++) {
					source[i] = ee[i];
				}
			}
			else if (j == num_x / 2) {
				for (int i = 0; i < numActiveOmega; i++) {
					source[i + numOmX - numActiveOmega] = ee[i + num_t * j];
				}
			}
			else {
				for (int i = 0; i < num_t; i++) {
					source[i + (j - 1) * num_t + numActiveOmega] = ee[i + num_t * j];
				}
			}
		}
	}
	else {
		for (int i = 0; i < num_t; i++) {
			double s = (-domain_t/2.0) + ht * ((double)i + 1);
			ee[i] = pparams->A0 * (sqrt(1.0 - pparams->relativeIntensity) * exp(-2.0 * log(2.0) * pow(s, 2) / (pow(pparams->pulseDuration, 2))) * cos(pparams->omega0 * s)
							+ sqrt(pparams->relativeIntensity) *       exp(-8.0 * log(2.0) * pow(s, 2) / (pow(pparams->pulseDuration, 2))) * cos(2.0 * pparams->omega0 * s + pparams->relativePhase));
			//ee[i] = 0.0;
		}
		writeInputEfield(ee);

		fftw_execute(e_f);
		normalizeFFT(ee, fftnorm);
		for (int i = 0; i < numActiveOmega; i++) {
			source[i] = ee[i];
		}
	}
	if(fetestexcept(FE_OVERFLOW | FE_INVALID | FE_DIVBYZERO)) raise(SIGFPE);
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
	//uniform_real_distribution<double> dis(-NOISE_MAGNITUDE, NOISE_MAGNITUDE);
	uniform_real_distribution<double> dis(0.0, NOISE_MAGNITUDE);
	
	/* complex<double>* integral1 = (complex<double>*)malloc(sizeof(complex<double>) * numActiveOmega);
	complex<double>* integral2 = (complex<double>*)malloc(sizeof(complex<double>) * numActiveOmega);
	complex<double>* Am_guess1 = (complex<double>*)malloc(sizeof(complex<double>) * numActiveOmega);
	complex<double>* Am_guess2 = (complex<double>*)malloc(sizeof(complex<double>) * numActiveOmega);

	for (int i = 0; i < numActiveOmega; i++)
	{
		integral1[i] = 0.0;
		integral2[i] = 0.0;
	}

	// Set the initial guess with y
	for (int k = 0; k < sizeRoot/2; k++){
		gsl_vector_set(u, k, odeObj->y[k + 2*numActiveOmega + freqLowerCutoff]);
		gsl_vector_set(u, k + sizeRoot/2, odeObj->y[k + 3*numActiveOmega + freqLowerCutoff]);
	}

	// Use integral condition to inform guess
    rootObj->setOutParam(1); // Turn output on (1) or off (0)
    rootObj->setIntCond(1); // Set whether to calculate integral
	rootObj->integral = integral1;
	mapG(u, rootObj, tmp);

	for (int k = 0; k < numActiveOmega; k++) {
		Am_guess1[k] = odeObj->y[k + 2 * numActiveOmega] + 1.0i * odeObj->y[k + 3 * numActiveOmega];
	}

	// Set the initial guess with y and uniform r.v.
	for (int k = 0; k < sizeRoot/2; k++){
		//gsl_vector_set(u, k, y[k + 2*numActiveOmega + freqLowerCutoff] + dis(gen) * INITIAL_GUESS_SEED_VALUE);
		//gsl_vector_set(u, k + sizeRoot/2, y[k + 3*numActiveOmega + freqLowerCutoff] + dis(gen) * INITIAL_GUESS_SEED_VALUE);
		gsl_vector_set(u, k, (1.0 + dis(gen)) * odeObj->y[k + 2*numActiveOmega + freqLowerCutoff]);
		gsl_vector_set(u, k + sizeRoot/2, (1.0 + dis(gen)) * odeObj->y[k + 3*numActiveOmega + freqLowerCutoff]);
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
		// -- Relative noise
		//gsl_vector_set(u, k, (1.0 + dis(gen)) * odeObj->y[k + 2*numActiveOmega + freqLowerCutoff]);
		//gsl_vector_set(u, k + sizeRoot/2, (1.0 + dis(gen)) * odeObj->y[k + 3*numActiveOmega + freqLowerCutoff]);
		//gsl_vector_set(u, k, (1.0 + dis(gen)) * gsl_vector_get(u, k));
		//gsl_vector_set(u, k, (1.0 + dis(gen)) * gsl_vector_get(u, k + sizeRoot/2));

		// -- Pure white noise
		//gsl_vector_set(u, k, dis(gen));
		//gsl_vector_set(u, k + sizeRoot/2, dis(gen));

		// -- Colored noise
		//gsl_vector_set(u, k, dis(gen) * odeObj->y[k + freqLowerCutoff]);
		//gsl_vector_set(u, k + sizeRoot/2, dis(gen) * odeObj->y[k + numActiveOmega + freqLowerCutoff]);
		//gsl_vector_set(u, k, dis(gen) * odeObj->y[k + 2*numActiveOmega + freqLowerCutoff]);
		//gsl_vector_set(u, k + sizeRoot/2, dis(gen) * odeObj->y[k + 3*numActiveOmega + freqLowerCutoff]);
		gsl_vector_set(u, k, dis(gen) * ((sizeRoot/2 - k) / (sizeRoot/2)) );
		gsl_vector_set(u, k + sizeRoot/2, dis(gen) * ((sizeRoot/2 - k) / (sizeRoot/2)));
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
