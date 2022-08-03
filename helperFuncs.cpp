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
	for (int its = 0; its < 5; its++) {
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
	int itSecant = 2;
	double f2norm;
	//ODEParams *odeObj = rootObj->getODEparams();
	gsl_vector *tmp = gsl_vector_alloc(sizeRoot);
	//complex<double> map1, map2; 
	double initGuess;

	// Create pseudo-random number generator for initial guess
	/* random_device rd;
	mt19937 gen(rd()); */
	default_random_engine gen;
	double rngLower, rngUpper;
	rngLower = -NOISE_MAGNITUDE;
	//rngLower = 0.0;
	rngUpper = NOISE_MAGNITUDE;
	uniform_real_distribution<double> dis(rngLower, rngUpper);
	printf(" Noise is generated from uniform distribution on [%.1e,%.1e)\n", rngLower, rngUpper);
	printf(" Noise type : Colored\n\n");
	
	complex<double>* GA1 = (complex<double>*)malloc(sizeof(complex<double>) * numActiveOmega);
	complex<double>* GA2 = (complex<double>*)malloc(sizeof(complex<double>) * numActiveOmega);
	complex<double>* Am_guess1 = (complex<double>*)malloc(sizeof(complex<double>) * numActiveOmega);
	complex<double>* Am_guess2 = (complex<double>*)malloc(sizeof(complex<double>) * numActiveOmega);

	for (int i = 0; i < numActiveOmega; i++)
	{
		GA1[i] = 0.0;
		GA2[i] = 0.0;
	}

	// Set the initial guess with y
	for (int k = 0; k < sizeRoot/2; k++){
		gsl_vector_set(u, k, odeObj->y[k + 2*numActiveOmega + freqLowerCutoff]);
		gsl_vector_set(u, k + sizeRoot/2, odeObj->y[k + 3*numActiveOmega + freqLowerCutoff]);
	}

	for (int k = 0; k < sizeRoot/2; k++) {
		//Am_guess1[k] = odeObj->y[k + 2 * numActiveOmega] + 1.0i * odeObj->y[k + 3 * numActiveOmega];
		Am_guess1[k] = gsl_vector_get(u, k) + 1.0i * gsl_vector_get(u, k + sizeRoot / 2);
	}

	// Use integral condition to inform guess
    rootObj->setOutParam(1); // Turn output on (1) or off (0)
	rootObj->getODEparams()->printAtomicProfile = 1;
    //rootObj->setIntCond(1); // Set whether to calculate integral
	//rootObj->integral = GA1;
	f2norm = mapG(u, rootObj);
	rootObj->getODEparams()->printAtomicProfile = 0;
	printf(" Iteration =  1, ||G||_2 = %.7e\n", f2norm);

	for (int k = 0; k < numActiveOmega; k++) {
		GA1[k] = odeObj->y[k + 2 * numActiveOmega] + 1.0i * odeObj->y[k + 3 * numActiveOmega];
	}

	// Set the initial guess with y and uniform r.v.
	for (int k = 0; k < sizeRoot/2; k++){
		//gsl_vector_set(u, k, y[k + 2*numActiveOmega + freqLowerCutoff] + dis(gen) * INITIAL_GUESS_SEED_VALUE);
		//gsl_vector_set(u, k + sizeRoot/2, y[k + 3*numActiveOmega + freqLowerCutoff] + dis(gen) * INITIAL_GUESS_SEED_VALUE);
		
		// Noisy linear solution
		gsl_vector_set(u, k, (1.0 + dis(gen)) * gsl_vector_get(u, k));
		gsl_vector_set(u, k + sizeRoot/2, (1.0 + dis(gen)) * gsl_vector_get(u, k + sizeRoot/2));

		// Colored noise
		//gsl_vector_set(u, k, dis(gen) * ((sizeRoot/2 - k) / (sizeRoot/2)) );
		//gsl_vector_set(u, k + sizeRoot/2, dis(gen) * ((sizeRoot/2 - k) / (sizeRoot/2)));

	}

	for (int k = 0; k < sizeRoot/2; k++) {
		//Am_guess2[k] = odeObj->y[k + 2 * numActiveOmega] + 1.0i * odeObj->y[k + 3 * numActiveOmega];
		Am_guess2[k] = gsl_vector_get(u, k) + 1.0i * gsl_vector_get(u, k + sizeRoot / 2);
	}

	// Get second guess for secant method
    rootObj->setItNum(2);
    //rootObj->setIntCond(2);
	//rootObj->integral = GA2;
	f2norm = mapG(u, rootObj);
	printf(" Iteration =  2, ||G||_2 = %.7e\n", f2norm);

	for (int k = 0; k < numActiveOmega; k++) {
		GA2[k] = odeObj->y[k + 2 * numActiveOmega] + 1.0i * odeObj->y[k + 3 * numActiveOmega];
	}

	// Do Secant update
	

	//cout << "Finished iteration 3" << endl;
	//if(fetestexcept(FE_OVERFLOW | FE_INVALID | FE_DIVBYZERO)) raise(SIGFPE); 

	for (int it = 3; it <= itSecant; it++) {

		for (int k = 0; k < sizeRoot/2; k++){
			//map1 = Am_guess1[k + freqLowerCutoff] + integral1[k + freqLowerCutoff];
			//map2 = Am_guess2[k + freqLowerCutoff] + integral2[k + freqLowerCutoff];
			if (real(GA2[k] - GA1[k]) < DBL_MIN) {
				gsl_vector_set(u, k, real(Am_guess1[k + freqLowerCutoff]));
			}
			else {
				//initGuess = (Am_guess2[k + freqLowerCutoff] * GA1[k] - Am_guess1[k + freqLowerCutoff] * GA2[k]) 
				// / (GA1[k] - GA2[k]);
				initGuess = (real(Am_guess2[k + freqLowerCutoff]) * real(GA1[k]) - real(Am_guess1[k + freqLowerCutoff]) * real(GA2[k])) 
				 / real(GA1[k] - GA2[k]);
				gsl_vector_set(u, k, initGuess);
				
			}

			if (imag(GA2[k] - GA1[k]) < DBL_MIN) {
				gsl_vector_set(u, k + sizeRoot/2, imag(Am_guess1[k + freqLowerCutoff]));
			}
			else {
				initGuess = (imag(Am_guess2[k + freqLowerCutoff]) * imag(GA1[k]) - imag(Am_guess1[k + freqLowerCutoff]) * imag(GA2[k])) 
				 / imag(GA1[k] - GA2[k]);
				gsl_vector_set(u, k + sizeRoot/2, initGuess);
			}
			
		}

		// Move latest update to previous guess
		for (int k = 0; k < numActiveOmega; k++) {
			Am_guess2[k] = Am_guess1[k]; 
			GA2[k] = GA1[k];
			//integral2[k] = integral1[k];
			//integral1[k] = 0.0;
		}

		for (int k = 0; k < sizeRoot/2; k++){
			Am_guess1[k] = gsl_vector_get(u, k) + 1.0i * gsl_vector_get(u, k + sizeRoot / 2);
		}
		
		/* for (int k = 0; k < sizeRoot/2; k++) {
			//Am_guess2[k] = odeObj->y[k + 2 * numActiveOmega] + 1.0i * odeObj->y[k + 3 * numActiveOmega];
			Am_guess2[k] = gsl_vector_get(u, k) + 1.0i * gsl_vector_get(u, k + sizeRoot / 2);
		} */

		// Do another iteration
		rootObj->setItNum(it);
		//rootObj->setIntCond(it); // Set whether to calculate integral
		//rootObj->integral = GA1;
		f2norm = mapG(u, rootObj);

		printf(" Iteration = %2d, ||G||_2 = %.7e\n", it, f2norm);

		for (int k = 0; k < sizeRoot/2; k++){
			GA1[k] = odeObj->y[k + 2 * numActiveOmega] + 1.0i * odeObj->y[k + 3 * numActiveOmega];
		}
	} 

	// Free some memory
	free(Am_guess1);
	free(Am_guess2);
	//free(integral1);
	//free(integral2);
	free(GA1);
	free(GA2);
	gsl_vector_free(tmp);
	
	
	/* cout << "Finished iteration 4" << endl;
	if(fetestexcept(FE_OVERFLOW | FE_INVALID | FE_DIVBYZERO)) raise(SIGFPE); */
	/* for (int k = 0; k < sizeRoot/2; k++){
		// -- Relative noise
		//gsl_vector_set(u, k, (1.0 + dis(gen)) * odeObj->y[k + 2*numActiveOmega + freqLowerCutoff]);
		//gsl_vector_set(u, k + sizeRoot/2, (1.0 + dis(gen)) * odeObj->y[k + 3*numActiveOmega + freqLowerCutoff]);
		//gsl_vector_set(u, k, (1.0 + dis(gen)) * gsl_vector_get(u, k));
		//gsl_vector_set(u, k + sizeRoot/2, (1.0 + dis(gen)) * gsl_vector_get(u, k + sizeRoot/2));

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
	} */


	rootObj->setItNum(itSecant + 1);
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
	//double epsrel = GSL_ROOT3_DBL_EPSILON;
	const double bU = pow(GSL_DBL_EPSILON, 1.0/4.0);
	const double bl = pow(GSL_DBL_EPSILON, 3.0/4.0);
	const double br = pow(GSL_DBL_EPSILON, 7.0/8.0);

	const double G0 = rootObj->getGnorm();
	complex<double> *linearRefl = rootObj->getLinSol();
	double *fac = rootObj->getDerivFactors();
	//gsl_vector *Am_f = gsl_vector_alloc(nRoot);
	//gsl_vector *Am_b = gsl_vector_alloc(nRoot);
	//gsl_vector_memcpy(Am_f, Am_in);
	//gsl_vector_memcpy(Am_b, Am_in);

	#pragma omp parallel
	{	
		//printf("Creating local GSL objects in gradient..\n");
		ODEParams private_odeObj = *(rootObj->getODEparams());
    	RootParams private_rootObj(&private_odeObj, nRoot);
		private_rootObj.setItNum(rootObj->getItNum());
		private_rootObj.setLinSol(linearRefl);

		int num_threads = omp_get_num_threads();
		int col_frac = nRoot / num_threads;
		int thread_id = omp_get_thread_num();
		int last_index;

		gsl_vector *Am_f;

		if (thread_id == num_threads - 1)
			last_index = nRoot;
		else
			last_index = (thread_id + 1) * col_frac;

		#pragma omp critical
		{
			//printf("Allocating copy of guess in gradient..\n");
			Am_f = gsl_vector_alloc(nRoot);
		}
		gsl_vector_memcpy(Am_f, Am_in);
		//printf("Copied guess in gradient..\n");

		//for (int j = 0; j < nRoot; j++){
		for (int j = thread_id * col_frac; j < last_index; j++) {
			double s, diff2, dGij, deltaA;
			double Aj = gsl_vector_get(Am_in, j);
			deltaA = fac[j] * abs(Aj);
			//double deltaA = epsrel;

			if (deltaA == 0) {
				deltaA = GSL_SQRT_DBL_EPSILON;
			}

			gsl_vector_set(Am_f, j, Aj + deltaA);
			//gsl_vector_set(Am_b, j, Aj - deltaA);

			double Gf = mapG(Am_f, &private_rootObj);
			//double Gb = mapG(Am_b, rootparams);
			double diff = Gf-G0;
			dGij = diff / deltaA;
			
			// Set scale, s, as max of two evaluations
			if (Gf > G0) s = Gf;
			else s = G0;
			
			// Adjust step size for future steps
			if (abs(diff) > s * br) {
				fac[j] = fac[j] * GSL_SQRT_DBL_EPSILON;
			}
			if (abs(diff) > s * br && abs(diff) < s * bl) {
				fac[j] = fac[j] / GSL_SQRT_DBL_EPSILON;
			}
			// If difference incorrect scale, then recompute the column
			if (diff < s * br) {
				deltaA = sqrt(fac[j]) * abs(Aj);
				Gf = mapG(Am_f, &private_rootObj);

				diff2 = Gf - G0;

				// Set the scale with new evaluation
				if (Gf > G0) s = Gf;
				else s = G0;

				// If step is acceptable, then accept new step and adjust factor
				if (sqrt(fac[j]) * abs(diff2) > abs(diff)) {
					dGij = diff2 / deltaA;

					if (abs(diff) > s * br) {
						fac[j] = fac[j] * GSL_SQRT_DBL_EPSILON;
					}
					if (abs(diff2) > s * br && abs(diff2) < s * bl) {
						fac[j] = fac[j] / GSL_SQRT_DBL_EPSILON;
					}
				}
			}

			gsl_vector_set(dG, j, dGij);
			//gsl_vector_set(dG, j, (Gf-Gb)/(2.0*deltaA));

			gsl_vector_set(Am_f, j, Aj);
			//gsl_vector_set(Am_b, j, Aj);
		}
	}
}

void GdG(const gsl_vector *Am_in, void *rootparams, double *G, gsl_vector *dG) {
	*G = mapG(Am_in, rootparams);
	dGdA(Am_in, rootparams, dG);
}
