#include "Simulation.h"

double min_mapG(const gsl_vector *ym_guess, void *simPtr);

void dGdA(const gsl_vector *Am_in, void *rootparams, gsl_vector *dG) {
	RootParams *rootObj = reinterpret_cast<RootParams*>(rootparams);
	const int nRoot = rootObj->getSizeRoot();
	//double epsrel = 1e-4;
	double epsrel = GSL_SQRT_DBL_EPSILON; // 1.49e-8
	//double epsrel = GSL_ROOT3_DBL_EPSILON;

	const double G0 = rootObj->getGnorm();
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
			double Aj = gsl_vector_get(Am_in, j);
			double deltaA = epsrel * abs(Aj);
			//double deltaA = epsrel;

			if (deltaA == 0) {
				deltaA = epsrel;
			}

			gsl_vector_set(Am_f, j, Aj + deltaA);
			//gsl_vector_set(Am_b, j, Aj - deltaA);

			double Gf = mapG(Am_f, &private_rootObj);
			//double Gb = mapG(Am_b, rootparams);

			gsl_vector_set(dG, j, (Gf-G0)/(deltaA));
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




void Simulation::fdfMin_iterateBPPE() {
	// Find the size of the problem with omega cutoffs
	int sizeRoot;
	ODEParams myODEParams(m_numT, m_omegaArray);
	if (numDimensionsMinusOne == 1) {
		sizeRoot = 2 * m_numX * (m_freqUpperCutoff - m_freqLowerCutoff + 1);
	}
	else {
		sizeRoot = 2 * (freqUpperCutoff - freqLowerCutoff + 1);
	}

    // Create the gslParams objects
    RootParams myRootParams(&myODEParams, sizeRoot);

    // Fill material-specific parameters with values from first layer
	myODEParams.fillParams(m_structure->m_layers.begin()->getMaterial());

	// Initialize time and status variables
	int status; 
	double dxnorm, xnorm, fnorm;
	double fval, simplexSize, gradnorm;
	double nonlinear_time_initial, nonlinear_time, nonlinear_time_tmp, nonlinear_time_total;

	char timeLogFname[STRING_BUFFER_SIZE];
	snprintf(timeLogFname, sizeof(char) * STRING_BUFFER_SIZE, "%stimeLog.txt", SIM_DATA_OUTPUT);
	FILE *localLogFile = fopen(timeLogFname, "a");

	// Initialize objects for condition number checking
	/* gsl_matrix *U = gsl_matrix_alloc(sizeRoot, sizeRoot);
	gsl_matrix *V = gsl_matrix_alloc(sizeRoot, sizeRoot);
	gsl_vector *singularValues = gsl_vector_alloc(sizeRoot);
	gsl_vector *work = gsl_vector_alloc(sizeRoot); */

	// Initialize multiroot objects
	printf("Allocating multiroot solver\n");
	//const gsl_multimin_fminimizer_type *fminType = gsl_multimin_fminimizer_nmsimplex2;
    //gsl_multimin_fminimizer *gslSolver;

	const gsl_multimin_fdfminimizer_type *fdfminType = gsl_multimin_fdfminimizer_vector_bfgs2;
	//const gsl_multimin_fdfminimizer_type *fdfminType = gsl_multimin_fdfminimizer_conjugate_fr;
	gsl_multimin_fdfminimizer *gslSolver_fdf;

	nonlinear_time_initial = omp_get_wtime();
	//gslSolver = gsl_multimin_fminimizer_alloc(fminType, sizeRoot);
	gslSolver_fdf = gsl_multimin_fdfminimizer_alloc(fdfminType, sizeRoot);

	fprintf(localLogFile, "Time spent allocating multiroot solver : %f [s]\n", omp_get_wtime() - nonlinear_time_initial);
	nonlinear_time_initial = omp_get_wtime();

	
	// ------------- Initial Guess Finding -------------
	printf("Initializing array y..\n");
	myODEParams.initializeY(m_sourceLeft);
	printf("Doing linear problem..\n");
	Simulation::doLinearProblem(&myODEParams, m_sourceLeft, m_sourceRight, *m_structure);
	/* if(fetestexcept(FE_OVERFLOW | FE_INVALID | FE_DIVBYZERO)) raise(SIGFPE);
	feclearexcept(FE_ALL_EXCEPT); */

	gsl_vector *u = gsl_vector_alloc(sizeRoot);
	printf("Generating initial guess..\n");
	Simulation::generateGuess(u, &myRootParams, &myODEParams);
	/* if(fetestexcept(FE_OVERFLOW | FE_INVALID | FE_DIVBYZERO)) raise(SIGFPE);
	feclearexcept(FE_ALL_EXCEPT); */

	fprintf(localLogFile, "Time spent finding initial guess : %f [s]\n", omp_get_wtime() - nonlinear_time_initial);
	fclose(localLogFile);
	nonlinear_time_initial = omp_get_wtime();
	
	gsl_vector *tmp = gsl_vector_alloc(sizeRoot);
	myRootParams.setOutParam(1);
	mapG(u, &myRootParams);
	myRootParams.setItNum(myRootParams.getItNum() + 1);
	myRootParams.setOutParam(0);

	// ---------------------------------------------------

	// Tell GSL multiroot the function and initial guess
	printf("Setting multiroot function\n");
	fflush(stdout);
	//gsl_multimin_function minFunc;
	gsl_multimin_function_fdf minFunc;
	minFunc.n = sizeRoot;
	minFunc.f = &min_mapG;
	minFunc.params = &myRootParams;
	minFunc.df = &dGdA;
	minFunc.fdf = &GdG;
	nonlinear_time_tmp = omp_get_wtime();

	//gsl_vector *initStepSizes = gsl_vector_alloc(sizeRoot);
	//gsl_vector_set_all(initStepSizes, minInitStep);
	//status = gsl_multimin_fminimizer_set (gslSolver, &minFunc, u, initStepSizes);
	status = gsl_multimin_fdfminimizer_set (gslSolver_fdf, &minFunc, u, minInitStep, minTol);

	if (status == GSL_EBADFUNC) {
			printf("\nERROR: GSL fsolver returned GSL_EBADFUNC. The iteration scheme encountered a singular point.\n");
	}
	if (status == GSL_ESING) {
		printf("\nERROR: GSL fsolver returned GSL_ESING. The first estimate of the Jacobian is singular.\n");
	}	

	nonlinear_time = omp_get_wtime() - nonlinear_time_tmp;
	printf("Finished setting multiroot function in %.2f seconds.\n", nonlinear_time);

	localLogFile = fopen(timeLogFname, "a");
	fprintf(localLogFile, "Initializing fsolver took : %f [s]\n\n", omp_get_wtime() - nonlinear_time_initial);
	fprintf(localLogFile, "===============================================================================================\n");
	fprintf(localLogFile, "| Iteration Number | Time in seconds  | 2-norm of gradient  |   2-norm of dx   | 2-norm of map |\n");
	fprintf(localLogFile, "===============================================================================================\n");
	fclose(localLogFile);
	nonlinear_time_initial = omp_get_wtime();

	// Compute the condition number of the Jacobian
	/* gsl_matrix *Jac = (hybrid_state_t)(s->state)->J;
	gsl_matrix_memcpy(U, Jac);
	gsl_linalg_SV_decomp(U, V, singularValues, work); */

    // Comment to toggle output on every iteration
	myRootParams.setOutParam(0);
	do
	{
		printf("Starting iteration %d\n", myRootParams.getItNum());
		if (myRootParams.getItNum() % outputInterval == 0) {
			myRootParams.setOutParam(1);
		}
		else {
			myRootParams.setOutParam(0);
		}
		fflush(stdout);
		nonlinear_time_tmp = omp_get_wtime();
		//status = gsl_multimin_fminimizer_iterate(gslSolver);
		status = gsl_multimin_fdfminimizer_iterate(gslSolver_fdf);

		//if(fetestexcept(FE_ALL_EXCEPT & ~FE_INEXACT)) raise(SIGFPE);
		if (status == GSL_EBADFUNC) {
			printf("\nERROR: GSL fsolver returned GSL_EBADFUNC. The iteration scheme encountered a singular point.\n");
		}
		if (status == GSL_ENOPROG) {
			printf("\nERROR: GSL fsolver returned GSL_ENOPROG. The iteration scheme is not making progress.\n");
		}
		if (status == GSL_ENOPROGJ) {
			printf("\nERROR: GSL fsolver returned GSL_ENOPROGJ. The iteration scheme is not making progress.\n");
		}
		if (status) break;

		//simplexSize = gsl_multimin_fminimizer_size (gslSolver);
		//status = gsl_multimin_test_size (simplexSize, minStopCon);
		status = gsl_multimin_test_gradient (gslSolver_fdf->gradient, minStopCon);
		gradnorm = gsl_blas_dnrm2(gsl_multimin_fdfminimizer_gradient(gslSolver_fdf));
        dxnorm = gsl_blas_dnrm2(gsl_multimin_fdfminimizer_dx(gslSolver_fdf));

		nonlinear_time = omp_get_wtime() - nonlinear_time_tmp;
		printf("Iteration %d completed in %.3f seconds.\n", myRootParams.getItNum(), nonlinear_time);
		
		if (myRootParams.getItNum() % outputInterval == 0) {
			localLogFile = fopen(timeLogFname, "a");
			fprintf(localLogFile, "|%18d|%18.3f|%18.5e|%18.5e|%18.5e|\n", myRootParams.getItNum(), nonlinear_time, gradnorm, dxnorm, gsl_multimin_fdfminimizer_minimum(gslSolver_fdf));
			//fprintf(localLogFile, "|%18d|%18.3f|%18.5e|%18.5e|\n", myRootParams.getItNum(), nonlinear_time, simplexSize, gsl_multimin_fminimizer_minimum(gslSolver));
			fclose(localLogFile);
		}	

		myRootParams.setItNum(myRootParams.getItNum() + 1);
		//fflush(stdout);
	}
	while (status == GSL_CONTINUE && myRootParams.getItNum() < maxIter);

	nonlinear_time_total = omp_get_wtime() - nonlinear_time_initial;
	printf("==========================================================\n");
	printf("Multiroot solver completed in %.2f seconds.\n\n", nonlinear_time_total);

	localLogFile = fopen(timeLogFname, "a");
	fprintf(localLogFile, "============================================================================\n");
	fprintf(localLogFile, "Quasi-Newton scheme has stopped after : %f [s]\n", nonlinear_time_total);
	fprintf(localLogFile, "Scheme return code : %d \n", status);
	fclose(localLogFile);


	gradnorm = gsl_blas_dnrm2(gsl_multimin_fdfminimizer_gradient(gslSolver_fdf));
    dxnorm = gsl_blas_dnrm2(gsl_multimin_fdfminimizer_dx(gslSolver_fdf));
	//printf("  simplex size = %.7e\n", simplexSize);
	printf("  gradient norm = %.7e\n", gradnorm);
    printf("        dx norm = %.7e\n", dxnorm);
    //printf("  f min = %.7e\n\n", gsl_multimin_fminimizer_minimum(gslSolver));
    printf("          f min = %.7e\n\n", gsl_multimin_fdfminimizer_minimum(gslSolver_fdf));


	// Run the map one last time to output spectra
	printf("Performing final iteration [%d] with output enabled..\n", myRootParams.getItNum());
	myRootParams.setOutParam(1);
	//mapG(gslSolver->x, &myRootParams);
	mapG(gslSolver_fdf->x, &myRootParams);

	// Free solver memory
	printf("Freeing solver memory.\n");
	gsl_vector_free(u);
    //gsl_multiroot_fsolver_free(s);
	//gsl_multimin_fminimizer_free(gslSolver);
	gsl_multimin_fdfminimizer_free(gslSolver_fdf);
	printf("Finished freeing solver memory.\n");
	
}


void Simulation::readParamfile(char *inFile) {
	char *paramFileBuffer = readParmetersFileToBuffer(inFile);

	// Load in domain parameters
	m_numT = getIntParameterValueByName("numTimePoints");
	m_domT = getDoubleParameterValueByName("timeDomainSize");
	m_numX = getIntParameterValueByName("numTransversePoints");
	m_domX = getDoubleParameterValueByName("transverseDomainSize");

	m_freqLowerCutoff = getIntParameterValueByName("omegLowerCutoff");
	m_freqUpperCutoff = getIntParameterValueByName("omegUpperCutoff");
	m_numActiveOmega = m_numT / 2 + 1;
	m_numOmX = m_numT*m_numX / 2 + 2;

	m_rho0 = getDoubleParameterValueByName("initialEDensity");

	// Read in optimization/root parameters
	m_maxIter = getIntParameterValueByName("multiminIterMax");
	m_outputInterval = getIntParameterValueByName("multiminOutputInterval");
	m_minInitStep = getDoubleParameterValueByName("multiminInitStep");
	m_minStopCon = getDoubleParameterValueByName("multiminStopCon");
	m_minTol = getDoubleParameterValueByName("multiminTol");
}
