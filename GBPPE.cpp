// bragg_linear.cpp : Defines the entry point for the console application.
// ACMS version 1

//#include "stdafx.h"

#include "BPPE.h"
#include "createLayers.h"
#include "helperFuncs.h"

using namespace std;

// Create GLOBAL structure and Material database
Structure myStructure("myFirstStructure");
Simulation mySimulation("andrewApplication1");
MaterialDB myMaterialsDB("myFirstMaterialDB");

// This block of vars were orignally inside main()
double sampleLayerThickness, I_0, A_0, tau, lambda_0, omega_0;
double LHSsourceLayerThickness, RHSbufferLayerThickness;
double twoColorSH_amplitude, twoColorSH_phase;
int num_t, freqLowerCutoff, freqUpperCutoff, numActiveOmega, numActiveOmega2, l_0, num_Threads;
double domain_t, zStepMaterial1, alpha_tukey;
vector<double> monitorZlocations;
int GSLerrorFlag, p, oFlag, VERBOSE;
double *omegaArray, *timeValuesArray, *kx;
complex<double> *sourceLeft, *sourceRight;
//fftw_plan nkForwardFFT, eFieldPlusForwardFFT, eFieldPlusBackwardFFT, eFieldMinusForwardFFT, eFieldMinusBackwardFFT, intBackwardFFT, npForwardFFT;
char *paramFileBuffer, SIM_DATA_OUTPUT[30];

int foundNaN = 0;
int delmeFLAG = 0;

int main(int argc, char *argv[])
{	
	// Start timer
	double dtime = omp_get_wtime();

	/* ============================================== */
	/* == Reading input file and processing input     */
	readGlobalParameters(argv[1]);

	pulseparam_type* sourceLeftParams = (pulseparam_type*)malloc(sizeof(pulseparam_type));
	pulseparam_type* sourceRightParams = (pulseparam_type*)malloc(sizeof(pulseparam_type));

	sourceLeftParams->A0 = sqrt(2.0 * I_0 / (epsilon_0*cLight));
	sourceLeftParams->omega0 = 2 * M_PI*cLight / lambda_0;
	sourceLeftParams->relativeIntensity = twoColorSH_amplitude;
	sourceLeftParams->relativePhase = twoColorSH_phase;
	sourceLeftParams->pulseDuration = tau;

	sourceRightParams->A0 = 0.0;
	sourceRightParams->omega0 = 2 * M_PI*cLight / lambda_0;
	sourceRightParams->relativeIntensity = twoColorSH_amplitude;
	sourceRightParams->relativePhase = twoColorSH_phase;
	sourceRightParams->pulseDuration = tau;
	
	/* ============================================== */
	/* == Do some setup and output simulation info    */

	// Set the signal handler for floating-point exceptions
	//feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);
	signal(SIGFPE, floatingPointExceptions); 
	signal(SIGNAN, floatingPointExceptions); 

	time_t now = time(0);
	char *datetime = ctime(&now);

	if (VERBOSE >= 0) {
		cout << "GBPPE code - ACMS Ver.1" << endl;
		cout << "The current date and time: " << datetime << endl;
		cout << "Verbosity = "<< VERBOSE << endl << endl;
		cout << "Working Directory = " << get_current_dir() << SIM_DATA_OUTPUT << endl << endl;
	}

	omp_set_num_threads(num_Threads);
	cout << "Num threads set to  = " << num_Threads << endl << "  Test: ";

	#pragma omp parallel 
	{
		#pragma omp critical
		{
			cout << omp_get_thread_num() << " ";
		}
	}
	cout << endl << endl;

	writeSimParameters();

	/* ============================================== */

	generateLayers(myMaterialsDB, myStructure);
    setupPointMonitorLocations(myMaterialsDB, myStructure);
	/// VERY Early termination
	//printf("!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!: VERY Early Termination\n"); exit(-1);

	printf("Structure generated and point monitor locations set..\n");
	
	sourceLeft = (complex<double>*)malloc(sizeof(complex<double>) * numActiveOmega);
	sourceRight = (complex<double>*)malloc(sizeof(complex<double>) * numActiveOmega);
	
	//y = (double*)malloc(sizeof(double) * 4 * numActiveOmega);

#ifndef _WIN64
	fftw_import_system_wisdom();
	char wisdomFilePath[STRING_BUFFER_SIZE];
	snprintf(wisdomFilePath, sizeof(char) * STRING_BUFFER_SIZE, "%smywisdomfile", SIM_DATA_OUTPUT);
	printf("\t\tREADING mywisdomfile: %s \n", wisdomFilePath);
    FILE* wisdomFile;
	if ((wisdomFile = fopen(wisdomFilePath, "r")) == NULL) { printf("WARNING: mywisdomfile not found\n"); }
	else {
		printf("mywisdomfile found importing wisdom...\n");
		fftw_import_wisdom_from_file(wisdomFile);
		fclose(wisdomFile);
	}
#endif
	
	omegaArray = (double*)malloc(sizeof(double)*num_t);
	kx = NULL;

	#ifndef _WIN64
	printf(" SAVING mywisdomfile...\n");
	if ((wisdomFile = fopen(wisdomFilePath, "w")) == NULL){
		printf("ERROR: CANNOT OPEN/ACCESS FILE - mywisdomfile\n");
		exit(-1);
	}
	else {
		fftw_export_wisdom_to_file(wisdomFile);
		fclose(wisdomFile);
	}
	#endif

	//FILE *fp;
	////errno_t err;
	//fp = fopen("n.dat", "w");


	fill_omg_k(omegaArray, kx, myMaterialsDB);
	DELME_ArgonDispersion(omegaArray);
	#ifdef DO_CONSTPLASMA
		DELME_AndrewPreformed(omegaArray, myMaterialsDB.getMaterialByName("PlasmaMat"));
	#endif
	#ifdef DO_ARGON_PLASMA
		DELME_AndrewPreformed(omegaArray, myMaterialsDB.getMaterialByName("Argon"));
	#endif
	createWindowFunc(alpha_tukey);
	//if (fp != NULL) { fclose(fp);  }
	

	/* ============================================== */
	/* == Generate sources */
	complex<double> *eFieldPlus = (complex<double>*)malloc(sizeof(complex<double>)*num_t);
	complex<double> *eFieldMinus = (complex<double>*)malloc(sizeof(complex<double>)*num_t);
	fftw_plan eFieldPlusForwardFFT = fftw_plan_dft_1d(num_t, reinterpret_cast<fftw_complex*>(&eFieldPlus[0]), reinterpret_cast<fftw_complex*>(&eFieldPlus[0]), FFTW_FORWARD, FFTW_WISDOM_TYPE );
	fftw_plan eFieldMinusForwardFFT = fftw_plan_dft_1d(num_t, reinterpret_cast<fftw_complex*>(&eFieldMinus[0]), reinterpret_cast<fftw_complex*>(&eFieldMinus[0]), FFTW_FORWARD, FFTW_WISDOM_TYPE );

	if (VERBOSE >=3) printf("Generating the right-hand side source.\n");
	generateTwoColorPulse(eFieldMinus, eFieldMinusForwardFFT, sourceRight, sourceRightParams);
	for (int i = 0; i < numActiveOmega; i++) {
		sourceRight[i] = 0.0;
	} 
	if (VERBOSE >=3) printf("Generating the left-hand side source.\n");
    generateTwoColorPulse(eFieldPlus, eFieldPlusForwardFFT, sourceLeft, sourceLeftParams);

	// Free time-domain fields and fft plans
	free(eFieldPlus);
	free(eFieldMinus);
	fftw_destroy_plan(eFieldPlusForwardFFT);
	fftw_destroy_plan(eFieldMinusForwardFFT);

	writeInputSpectrum(sourceLeft);
	/* ============================================== */
	
	std::string reldatpath = SIM_DATA_OUTPUT;
	myStructure.writeStructureLayoutToASCIIFile(reldatpath + "StructureLayout.txt");
	myStructure.writeStructureToDATFile(reldatpath + "Structure.dat");
	myMaterialsDB.writeMaterialDBToASCIIFile(reldatpath + "MaterialDatabase.txt");
	myStructure.writeBoundaryLayoutToASCIIFile(reldatpath + "BoundaryLayout.txt");

	/// Early termination
	//printf("############ WARNING ###########:  Early Termination\n"); exit(-1);
    printf("\n\n############ INFO ###########:  ENTERING THE NONLINEAR PART ##############################\n");

	char timeLogFname[STRING_BUFFER_SIZE];
	snprintf(timeLogFname, sizeof(char) * STRING_BUFFER_SIZE, "%stimeLog.txt", SIM_DATA_OUTPUT);
	FILE *tLogFile = fopen(timeLogFname, "w");
	fprintf(tLogFile, "Time spent before entering iterateBPPE() : %f [s]\n\n", omp_get_wtime() - dtime);
	fprintf(tLogFile, "============================================================================\n");
	fclose(tLogFile);

	/* ============================================== */
	/* == Primary computation */
	iterateBPPE();

	dtime = omp_get_wtime() - dtime;
	/* ============================================== */

	if(VERBOSE >= 0) { cout << "Time in seconds is " << dtime << endl; }

	tLogFile = fopen(timeLogFname, "a");
	fprintf(tLogFile, "\n==========================================================\n\n");
	fprintf(tLogFile, "Total time spent in program: %f [s]\n\n", dtime);
	fclose(tLogFile);
	//cout << "Num threads set to  = " << omp_get_num_threads() << endl << endl;
	cout << "The time of various steps have been recorded in the following file: " << timeLogFname << endl;
	
	cout << "Cleaning up memory further.." << endl;
	free(sourceLeft);
	free(sourceRight);
	

	cout << endl << "Exiting program.." << endl << endl;

    return 0;
}

void setupPointMonitorLocations(MaterialDB& theMaterialDB, Structure& theStructure)
{

	int num10ums = (myStructure.getThickness() - RHSbufferLayerThickness) / 10e-6;
	for (int n = 1; n < num10ums; n++){
		monitorZlocations.push_back(LHSsourceLayerThickness + n*10e-6);
	}
	monitorZlocations.push_back(myStructure.getThickness() - RHSbufferLayerThickness);

}



int mapG(const gsl_vector *ym_guess, void *rootparams, gsl_vector *f) {
    //rootparam_type *rparams = reinterpret_cast<rootparam_type*>(rootparams);
    RootParams *rootObj = reinterpret_cast<RootParams*>(rootparams);

	// Initialize GSL ODE objects
	const gsl_odeiv2_step_type * stepType = gsl_odeiv2_step_rkf45;
	gsl_odeiv2_step * gslStep = gsl_odeiv2_step_alloc(stepType, 4 * numActiveOmega);
	gsl_odeiv2_control * gslControl = gsl_odeiv2_control_y_new(ode_epsabs, ode_epsrel);
	gsl_odeiv2_evolve * gslEvolve = gsl_odeiv2_evolve_alloc(4 * numActiveOmega);

	gsl_odeiv2_system sys = { dAdz, NULL, (size_t)(4 * numActiveOmega), rootObj->getODEparams()};
	gsl_odeiv2_evolve_reset(gslEvolve);

	double *yloc = (rootObj->getODEparams())->y;

    double nonlinear_time_initial, nonlinear_time;

	int numzReports, numZsteps = 0;
	double zRight, zStepSize;
	double zPosition = 0.0;
	vector<double> endPoints;

    nonlinear_time_initial = omp_get_wtime();

	// Reset the left source and update guess
	updateGuess(yloc, sourceLeft, ym_guess, rootObj);
	
	if (rootObj->getOutParam() == 1) {
		write_out_eFieldAndSpectrumAtZlocation(rootObj->getItNum(), 
            0, yloc, 0.0, rootObj->getODEparams()->ee_m, 
            myStructure.m_layers.front().getMaterial().getK(), 
            rootObj->getODEparams()->em_b);
	}
    //cout << "  Going FORWARD through layers" << endl;
    for (std::list<Layer>::iterator lit = myStructure.m_layers.begin(); lit != myStructure.m_layers.end(); ++lit) {
        // Skip the LHS layer and the RHS layers
        if (lit->getLowSideBoundary() != NULL && lit->getHiSideBoundary() != NULL)
        {
            boundary(lit->getLowSideBoundary()->m_zPos, lit->getLowSideBoundary()->lowSideLayer->getMaterial().getK(), lit->getLowSideBoundary()->hiSideLayer->getMaterial().getK(), yloc);
            /* rootObj->getODEparams()->k = lit->getMaterial().getK();
            rootObj->getODEparams()->chi_2 = lit->getMaterial().getChi2();
            rootObj->getODEparams()->chi_3 = lit->getMaterial().getChi3();
            rootObj->getODEparams()->doPlasmaCalc = lit->getMaterial().getdoPlasmaCalc();
            rootObj->getODEparams()->mpi_k = lit->getMaterial().getmpi_k();
            rootObj->getODEparams()->mpi_sigmaK = lit->getMaterial().getmpi_sigmaK();
            rootObj->getODEparams()->ionE = lit->getMaterial().getIonizationEnergy(); */

			rootObj->getODEparams()->fillParams(lit->getMaterial());
			gsl_odeiv2_evolve_reset(gslEvolve);

            zStepSize = lit->getStepSize();
            zPosition = lit->getStartZpos();
            //zRight = lit->getEndZpos();

            for (std::size_t atZpos = 0; atZpos < monitorZlocations.size(); ++atZpos) {
                if (monitorZlocations[atZpos] > lit->getStartZpos() && monitorZlocations[atZpos] < lit->getEndZpos())
                {
                    endPoints.push_back(monitorZlocations[atZpos]);
                }
            }
            endPoints.push_back(lit->getEndZpos());

            numzReports = 0; 
            //numZsteps = 0;

            //if (VERBOSE >= 6) { cout << endl << " Doing Layer #" << lit->getlayerIDnum() << endl; }
            for (std::size_t atZpos = 0; atZpos < endPoints.size(); ++atZpos)  {
                zRight = endPoints[atZpos];
                while (zPosition < zRight) 
                {

                    GSLerrorFlag = gsl_odeiv2_evolve_apply(gslEvolve, gslControl, gslStep, &sys, &zPosition, zRight, &zStepSize, yloc);

					for (int k = 0; k < numActiveOmega; k++){
						if(yloc[k] != yloc[k]) foundNaN++;
						if(yloc[k + numActiveOmega] != yloc[k + numActiveOmega]) foundNaN++;
						if(yloc[k + 2 * numActiveOmega] != yloc[k + 2 * numActiveOmega]) foundNaN++;
						if(yloc[k + 3 * numActiveOmega] != yloc[k + 3 * numActiveOmega]) foundNaN++;
					}

					if (foundNaN >= 1) {
						printf("WARNING: found NaN\n");
						printf("   Info: While taking ode step at z = %.5e\n", zPosition);
						raise(SIGABRT);
					}
                    
					if (GSLerrorFlag == GSL_SUCCESS) {
                        numZsteps++;
                    }
                    else {
                        printf("error: driver returned %d\n", GSLerrorFlag);
                        break;
                    }

                    nonlinear_time = omp_get_wtime() - nonlinear_time_initial;
                    if ((int)(nonlinear_time / 20) > numzReports && rootObj->getOutParam() == 1) {
                        printf("  I = %d, step = %d, z = %.8g, t = %d s\n", rootObj->getItNum(), numZsteps, zPosition, (int)nonlinear_time);
                        numzReports++;
                    }

					if (rootObj->getIntCond() != 0) {
						integrate(zPosition, zStepSize, rootObj->getODEparams(), yloc, rootObj->integral);
					}
                }

				if (rootObj->getOutParam() == 1) {
					//printf("  Outputting Point Monitor file at Z location %d[nm]... \n", (int)round(zPosition * 1.0e9));
					write_multicolumnMonitor(rootObj->getItNum(), zPosition, yloc, rootObj->getODEparams());

				}
				//if(fetestexcept(FE_ALL_EXCEPT & ~FE_INEXACT)) raise(SIGFPE);
                
            }

            endPoints.clear();

            //integrate(zPosition, zStepSize, params, y, integral);
            
        }
        //  Finally do the lowside boundary of the last assumed Vacuum layer
        if (lit->getHiSideBoundary() == NULL) {
            boundary(lit->getLowSideBoundary()->m_zPos, lit->getLowSideBoundary()->lowSideLayer->getMaterial().getK(), lit->getLowSideBoundary()->hiSideLayer->getMaterial().getK(), yloc);
        }
    }

	for (int k = 0; k < rootObj->getSizeRoot()/2; k++){
		gsl_vector_set(f, k, yloc[k + 2 * numActiveOmega + freqLowerCutoff] - real(sourceRight[k + freqLowerCutoff]));
		gsl_vector_set(f, k + rootObj->getSizeRoot()/2, yloc[k + 3 * numActiveOmega + freqLowerCutoff] - imag(sourceRight[k + freqLowerCutoff]));
	}


	if (rootObj->getOutParam() == 1) {
		write_out_eFieldAndSpectrumAtZlocation(rootObj->getItNum(), 
        1, yloc, myStructure.getThickness(), rootObj->getODEparams()->ee_p, 
        myStructure.m_layers.back().getMaterial().getK(), 
        rootObj->getODEparams()->ep_b);
	}
    
	nonlinear_time = omp_get_wtime() - nonlinear_time_initial;
	//if (rootObj->output == 1) {
		//printf("Iteration %d completed in %.2f seconds with %d steps.\n", rootObj->itnum, nonlinear_time, numZsteps);
		//cout << "Iteration " << rootObj->itnum <<  " completed in " <<  nonlinear_time << "seconds with" << numZsteps << "steps." << endl;
		//fflush(stdout);
	//}

	for (int k = 0; k < numActiveOmega; k++){
		yloc[k + 2*numActiveOmega] = real(sourceRight[k]);
		yloc[k + 3*numActiveOmega] = imag(sourceRight[k]);
	}

	myStructure.doBackwardPassThroughAllBoundaries(yloc);

	// Freeing ODE memory
	gsl_odeiv2_control_free(gslControl);
    gsl_odeiv2_evolve_free(gslEvolve);
    gsl_odeiv2_step_free(gslStep);

	for (int k = 0; k < numActiveOmega; k++){
		if(yloc[k] != yloc[k]) foundNaN++;
		if(yloc[k + numActiveOmega] != yloc[k + numActiveOmega]) foundNaN++;
		if(yloc[k + 2 * numActiveOmega] != yloc[k + 2 * numActiveOmega]) foundNaN++;
		if(yloc[k + 3 * numActiveOmega] != yloc[k + 3 * numActiveOmega]) foundNaN++;
	}

	if (foundNaN >= 1) 
		return GSL_EBADFUNC;
	else
    	return GSL_SUCCESS;
}

void iterateBPPE()
{
	// Find the size of the problem with omega cutoffs
	int sizeRoot = 2*(freqUpperCutoff - freqLowerCutoff + 1);

    // Create the gslParams objects
    ODEParams myODEParams(num_t, omegaArray);
    RootParams myRootParams(&myODEParams, sizeRoot);
	//RootParams myRootParams(num_t, omegaArray, sizeRoot);

    // Fill material-specific parameters with values from first layer
    //myRootParams.getODEparams()->fillParams(myStructure.m_layers.begin()->getMaterial());
	myODEParams.fillParams(myStructure.m_layers.begin()->getMaterial());

	// Initialize time and status variables
	int status; 
	double dxnorm, xnorm, fnorm;
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
	const gsl_multiroot_fsolver_type *T;
    gsl_multiroot_fsolver *s;
	T = gsl_multiroot_fsolver_hybrid;

	nonlinear_time_initial = omp_get_wtime();
	s = gsl_multiroot_fsolver_alloc(T, sizeRoot);

	fprintf(localLogFile, "Time spent allocating multiroot solver : %f [s]\n", omp_get_wtime() - nonlinear_time_initial);
	nonlinear_time_initial = omp_get_wtime();

	
	// ------------- Initial Guess Finding -------------
	printf("Initializing array y..\n");
	//initializeY(myRootParams.getODEparams()->y, sourceLeft);
	myODEParams.initializeY(sourceLeft);
	printf("Doing linear problem..\n");
	doLinearProblem(&myODEParams, sourceLeft, sourceRight, myStructure);
	/* if(fetestexcept(FE_OVERFLOW | FE_INVALID | FE_DIVBYZERO)) raise(SIGFPE);
	feclearexcept(FE_ALL_EXCEPT); */

	gsl_vector *u = gsl_vector_alloc(sizeRoot);
	printf("Generating initial guess..\n");
	generateGuess(u, &myRootParams, &myODEParams);
	/* if(fetestexcept(FE_OVERFLOW | FE_INVALID | FE_DIVBYZERO)) raise(SIGFPE);
	feclearexcept(FE_ALL_EXCEPT); */

	fprintf(localLogFile, "Time spent finding initial guess : %f [s]\n", omp_get_wtime() - nonlinear_time_initial);
	nonlinear_time_initial = omp_get_wtime();

	gsl_vector *tmp = gsl_vector_alloc(sizeRoot);
	myRootParams.setOutParam(1);
	mapG(u, &myRootParams, tmp);
	myRootParams.setItNum(myRootParams.getItNum() + 1);
	myRootParams.setOutParam(0);

	// ---------------------------------------------------

	// Tell GSL multiroot the function and initial guess
	printf("Setting multiroot function\n");
	gsl_multiroot_function f = {&mapG, sizeRoot, &myRootParams};
	nonlinear_time_tmp = omp_get_wtime();
	gsl_multiroot_fsolver_set(s, &f, u);
	nonlinear_time = omp_get_wtime() - nonlinear_time_tmp;
	printf("Finished setting multiroot function in %.2f seconds.\n", nonlinear_time);

	fprintf(localLogFile, "Initializing fsolver took : %f [s]\n\n", omp_get_wtime() - nonlinear_time_initial);
	fprintf(localLogFile, "============================================================================\n");
	fprintf(localLogFile, "| Iteration Number | Time in seconds  | 2-norm of step  | 2-norm of map    |\n");
	fprintf(localLogFile, "============================================================================\n");
	nonlinear_time_initial = omp_get_wtime();

	// Compute the condition number of the Jacobian
	/* gsl_matrix *Jac = (hybrid_state_t)(s->state)->J;
	gsl_matrix_memcpy(U, Jac);
	gsl_linalg_SV_decomp(U, V, singularValues, work); */

    // Comment to toggle output on every iteration
	myRootParams.setOutParam(1);
	do
	{
		printf("Starting iteration %d\n", myRootParams.getItNum());
		fflush(stdout);
		nonlinear_time_tmp = omp_get_wtime();
		status = gsl_multiroot_fsolver_iterate(s);

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

		//status = gsl_multiroot_test_residual(s->f, root_epsabs);
		status = gsl_multiroot_test_delta(s->dx, s->x, root_epsabs, root_epsrel);
		//dxnorm = gsl_blas_dasum(s->dx);
		dxnorm = gsl_blas_dnrm2(s->dx);
		fnorm = gsl_blas_dnrm2(s->f);

		nonlinear_time = omp_get_wtime() - nonlinear_time_tmp;
		printf("Iteration %d completed in %.3f seconds.\n", myRootParams.getItNum(), nonlinear_time);
		
		fprintf(localLogFile, "|%18d|%18.3f|%18.5e|%18.5e|\n", myRootParams.getItNum(), nonlinear_time, dxnorm, fnorm);

		myRootParams.setItNum(myRootParams.getItNum() + 1);
		fflush(stdout);
	}
	while (status == GSL_CONTINUE && myRootParams.getItNum() < 25000);

	nonlinear_time_total = omp_get_wtime() - nonlinear_time_initial;
	printf("==========================================================\n");
	printf("Multiroot solver completed in %.2f seconds.\n\n", nonlinear_time_total);

	fprintf(localLogFile, "============================================================================\n");
	fprintf(localLogFile, "Quasi-Newton scheme has stopped after : %f [s]\n", nonlinear_time_total);
	fclose(localLogFile);

	
	//dxnorm = gsl_blas_dasum(s->dx);
	//xnorm = gsl_blas_dasum(s->x);
	//fnorm = gsl_blas_dasum(s->f);
	dxnorm = gsl_blas_dnrm2(s->dx);
	xnorm = gsl_blas_dnrm2(s->x);
	fnorm = gsl_blas_dnrm2(s->f);
	printf("  dx norm = %.7e\n", dxnorm);
	printf("  x norm = %.7e\n", xnorm);
	printf("  f norm = %.7e\n\n", fnorm);


	// Run the map one last time to output spectra
	printf("Performing final iteration with output enabled..\n");
	myRootParams.setOutParam(1);
	mapG(s->x, &myRootParams, s->f);

	// Free solver memory
	printf("Freeing solver memory.\n");
	gsl_vector_free(u);
    gsl_multiroot_fsolver_free(s);
	printf("Finished freeing solver memory.\n");
	
}




#define ANDREW_PREFORMED_METHOD1
void DELME_AndrewPreformed(double* omg, Material* mat) {
	//Material *plasmaMat;
	//plasmaMat = myMaterialsDB.getMaterialByName("PlasmaMat");

	double preformedDensity = rho_0;
	double omega_plasma = sqrt(pow(charge_e, 2) * preformedDensity / (epsilon_0 * mass_e));
	printf("WARNING: Performing preformed plasma test with density %.1e.\n The plasma frequency is then: %.8e\n", preformedDensity, omega_plasma);

	char kFile[STRING_BUFFER_SIZE];
	snprintf(kFile, sizeof(char) * STRING_BUFFER_SIZE, "%skz_plasma.dat", SIM_DATA_OUTPUT);
	FILE* fp;
	fp = fopen(kFile, "w");
	if (fp != NULL)
	{
		fprintf(fp, "# omega [rads/s]\tRe[kz]\tIm[kz]\n");
	}
	else {
		printf("Failed to open file '%s'\n", kFile);

		exit(-1);
	}

	complex<double> n0;
	//plasmaMat->m_k[0] = 0.0;
	//fprintf(fp, "%.7g\t%.17g\t%.17g \n", 0.0, real(plasmaMat->m_k[0]), imag(plasmaMat->m_k[0]));
	mat->m_k[0] = 0.0;
	fprintf(fp, "%.7g\t%.17g\t%.17g \n", 0.0, real(mat->m_k[0]), imag(mat->m_k[0]));
	for (int i = 1; i < numActiveOmega; i++)
	{
		n0 = mat->m_k[i] * cLight / omg[i];
		//mat->m_k[i] = omg[i] / cLight * sqrt(complex<double>(1.0 - pow(omega_plasma / omg[i], 2)));
		mat->m_k[i] = omg[i] / cLight * sqrt(1.0 - pow(omega_plasma / omg[i], 2) / (1.0 + 1.0i / (omg[i] * tauCollision)));
		fprintf(fp, "%.7g\t%.17g\t%.17g \n", omg[i], real(mat->m_k[i]), imag(mat->m_k[i]));
	}
	if (fp != NULL) { fclose(fp); }
}

void DELME_ArgonDispersion(double* omg) {
	Material *ArgonMat;
	ArgonMat = myMaterialsDB.getMaterialByName("Argon");

	complex<double> n0;
	double lambda;
	double raylighFactor = 32.0 * pow(M_PI, 3) / (3.0 * pow(num_atoms,2)) * num_atoms;
	double rayleighScattering = 0.0;
	
	char dispersionFile[STRING_BUFFER_SIZE];
	snprintf(dispersionFile, sizeof(char) * STRING_BUFFER_SIZE, "%sn_Argon.dat", SIM_DATA_OUTPUT);
	FILE* fp;
	fp = fopen(dispersionFile, "w");
	if (fp != NULL)
	{
		fprintf(fp, "# lambda [um]\tRe[n]\tIm[n]\n");
	}
	else {
		printf("Failed to open file '%s'\n", dispersionFile);

		exit(-1);
	}

	ArgonMat->m_k[0] = 0.0;
	for (int i = 1; i < numActiveOmega; i++)
	{
		lambda = 2 * M_PI*cLight / omg[i] * 1e6;
		//n0 = 1+2.50141e-3/(91.012-pow(lambda,-2))+5.00283e-4/(87.892-pow(lambda,-2))+5.22343e-2/(214.02-pow(lambda,-2));
		//n0 = 1+6.432135E-5+2.8606021E-2/(144-pow(lambda,-2)); // From Peck and Fisher (1964), valid in range: 0.47-2.06 um

		if (lambda < 0.47e-6) {
			n0 = 1.0 + 2.50141e-3/(91.012-pow(lambda,-2)) + 5.00283e-4/(87.892-pow(lambda,-2)) + 5.22343e-2/(214.02-pow(lambda,-2)); // From Bideau-Mehu (1981), valid in range: 0.140-0.568 um
		}
		else {
			n0 = 1.0 + 6.7867e-5 + 3.0182943e-2/(144.0-pow(lambda,-2)); // From Peck and Fisher (1964), valid in range: 0.47-2.06 um (0 deg C)
		}

		// Fitted dispersion (in Mathematica)
		//  n0 = 1.00026 + 1.37647e-5 / lambda - 1.86046e-6 / pow(lambda,2) + 3.78465e-7 / pow(lambda,3);
		// NEW FIT (Using Cauchy equation)
		n0 = 1.00027 + 2.37091e-6 / pow(lambda, 2) + 1.19746e-9 / pow(lambda, 4);


		//n0 += 1.0i * 8.0 * pow(M_PI, 3) / (3 * num_atoms * pow(lambda,4)) * pow(pow(n0, 2) - 1.0, 2);
		rayleighScattering = raylighFactor / pow(lambda*1e-6, 4) * pow(real(n0)-1.0, 2);
		//n0 += 1.0i * rayleighScattering;

		//n0 = 1.00026; // COMMENT OUT! THIS USED FOR TESTING

		ArgonMat->m_k[i] = omg[i] * n0 / cLight;

		fprintf(fp, "%.7g\t%.17g\t%.17g \n", lambda, real(n0), imag(n0));
	}
	if (fp != NULL) { fclose(fp); }

}


void writeSimParameters()
{

	char parametersFilePathName[STRING_BUFFER_SIZE];
	snprintf(parametersFilePathName, sizeof(char) * STRING_BUFFER_SIZE, "%sSimParameters.dat", SIM_DATA_OUTPUT);

	FILE* fp;
	//errno_t err;
	if (VERBOSE >= 2) { printf("Writing Simulation Parameter file %s\n", parametersFilePathName); }
	fp = fopen(parametersFilePathName, "w");
	if (fp != NULL)
	{
		//fprintf(fp, "cLight         \t%.17g\t/*speed of light */\n", cLight);
		//fprintf(fp, "epsilon_0      \t%.17g\t/*permittivity of free space */\n", epsilon_0);
		//fprintf(fp, "Znaught        \t%.17g\t/*impedance of free space */\n", Znaught);
		//fprintf(fp, "charge_e       \t%.17g\t/*electron charge */\n", charge_e);
		//fprintf(fp, "mass_e         \t%.17g\t/*electron mass */\n", mass_e);
		fprintf(fp, "u_Argon        \t%.17g\t/*argon potential */\n", u_Argon);
		fprintf(fp, "u_Hydrogen     \t%.17g\t/*hydrogen potential */\n", u_Hydrogen);
		fprintf(fp, "E_a            \t%.17g\t/*QST parameters */\n", E_a);
		fprintf(fp, "mu_a           \t%.17g\t/*FILL */\n", mu_a);
		fprintf(fp, "I_0            \t%.8g\t/*initial peak ensity */\n", I_0);
		fprintf(fp, "twoColorSH_amplitude \t%.17g\t/*two-color pulse: 2nd harmonic with half duration of fundamental */\n", twoColorSH_amplitude);
		fprintf(fp, "twoColorSH_phase \t%.17g\t/*phase shift of 2nd harmonic */\n", twoColorSH_phase);
		fprintf(fp, "tau             \t%.8g\t/*pulse duration fwhm */\n", tau);
		fprintf(fp, "lambda_0        \t%.17g\t/*central wavelength */\n", lambda_0);
		fprintf(fp, "omega_0         \t%.17e\t/*central angular frequency */\n", omega_0);
		//fprintf(fp, "waist_x         \t%.17g\t/*pulse-waist fwhm */\n", waist_x);
		fprintf(fp, "num_t           \t%.17g\t/*number of time pos */\n", (double)num_t);
		fprintf(fp, "num_x           \t%.17g\t/*number of x pos */\n", (double)num_x);
		fprintf(fp, "domain_t        \t%.17g\t/*time domain */\n", domain_t);
		fprintf(fp, "domain_x        \t%.17g\t/*x domain */\n", domain_x);
		fprintf(fp, "num_atoms       \t%.17g\t/*number of atoms in gas */\n", num_atoms);
		fprintf(fp, "rho_0           \t%.17g\t/* initial electron density */\n", rho_0);
		fprintf(fp, "j_e0            \t%.17g\t/* initial current density */\n", j_e0);
		fprintf(fp, "freqUpperCutoff \t%d\t/* upper frequency cut-off */\n", freqUpperCutoff);
		fprintf(fp, "freqLowerCutoff \t%d\t/*lower frequency cut-off */\n", freqLowerCutoff);
		//fprintf(fp, "num_iterations  \t%d\t/*number of BPPE iterations*/ */\n", num_iterations);
		fprintf(fp, "shift           \t%.17g\t/*FILL */\n", shift);
		fprintf(fp, "numDimensionsMinusOne \t%.17g\t/*(1+1) dimension (0) or (2+1) dimension (1) */\n", (double)numDimensionsMinusOne);
		//fprintf(fp, "plasmaOnOff     \t%d\t/*plasma on (1) or off (0) */\n", plasmaOnOff);
		fprintf(fp, "l_0               \t%.17g\t/*FILL */\n", (double)l_0);
		fprintf(fp, "numActiveOmega    \t%.17g\t/*FILL */\n", (double)numActiveOmega);
		//fprintf(fp, "numActiveOmega2   \t%.17g\t/*FILL */\n", (double)numActiveOmega2);
		fprintf(fp, "omegaPlasmaDamping \t%.17g\t/*plasma damping */\n", omegaPlasmaDamping);
		fprintf(fp, "tauCollision       \t%.17g\t/*mean collision time */\n", tauCollision);
		//fprintf(fp, "lengthSample       \t%.17g\t/*length of slab */\n", lengthSample);
		fprintf(fp, "LHSsourceLayerThickness \t%.17g\t/*distance from laser source to slab */\n", LHSsourceLayerThickness);
		fprintf(fp, "RHSbufferLayerThickness \t%.17g\t/*distance from slab to receiver */\n", RHSbufferLayerThickness);
		//fprintf(fp, "zRightHandSideOfSample \t%.17g\t/*FILL */\n", zRightHandSideOfSample);
		//fprintf(fp, "Z_4                     \t%.17g\t/*FILL */\n", Z_4);
		fprintf(fp, "sampleLayerThickness   \t%.17g\t/*Thickness of structure*/\n", sampleLayerThickness);
		//fprintf(fp, "numLayersInSample      \t%.17g\t/*number of half-periods in slab */\n", (double)numLayersInSample);
		fprintf(fp, "zStepMaterial1         \t%.17g\t/*propagation aZstep in material 1 */\n", zStepMaterial1);
		//fprintf(fp, "zStepMaterial2         \t%.17g\t/*propagation aZstep in material 2 */\n", zStepMaterial2);
		//fprintf(fp, "numZstepsMaterial1     \t%.17g\t/*number of steps in material 1 */\n", (double)numZstepsMaterial1);
		//fprintf(fp, "numZstepsMaterial2     \t%.17g\t/*number of steps in material 2 */\n", (double)numZstepsMaterial2);
		fprintf(fp, "n0_Vacuum           \t%.17g\t/*central index in material 0 */\n", n0_Vacuum);
		fprintf(fp, "n0_Material1           \t%.17g\t/*central index in material 1 */\n", n0_Material1);
		fprintf(fp, "n0_Material2           \t%.17g\t/*central index in material 2 */\n", n0_Material2);
		fprintf(fp, "n0_Material3           \t%.17g\t/*central index in material 3 */\n", n0_Material3);
		fprintf(fp, "n2_Material1           \t%.17g\t/*nonlinear index in material 1 */\n", n2_Material1);
		fprintf(fp, "n2_Material2           \t%.17g\t/*nonlinear index in material 2 */\n", n2_Material2);
		fprintf(fp, "chi_2                  \t%.17g\t/*nonlinear chi_2 */\n", chi_2);
		fprintf(fp, "A_0                    \t%.17g\t/*pulse peak amplitude */\n", A_0);
		fprintf(fp, "Keldysh                \t%.17g\t/*keldysh parameter */\n", omega_0 * sqrt(2.0 * u_Argon*charge_e) / A_0);
		fprintf(fp, "omegaPlasma            \t%.17g\t/*fully-ionized plasma frequency */\n", omegaPlasma);
		fprintf(fp, "Sellmeir_chi_1_1       \t%.17g\t/*FILL */\n", Sellmeir_chi_1_1);
		fprintf(fp, "Sellmeir_chi_1_2       \t%.17g\t/*FILL */\n", Sellmeir_chi_1_2);
		fprintf(fp, "Sellmeir_chi_1_3       \t%.17g\t/*FILL */\n", Sellmeir_chi_1_3);
		fprintf(fp, "Sellmeir_omega_1       \t%.17g\t/*FILL */\n", Sellmeir_omega_1);
		fprintf(fp, "Sellmeir_omega_2       \t%.17g\t/*FILL */\n", Sellmeir_omega_2);
		fprintf(fp, "Sellmeir_omega_3       \t%.17g\t/*FILL */\n", Sellmeir_omega_3);

		// Print the GSL ODE parameters
		fprintf(fp, "ode_epsabs       \t%.17g\t/*absolute error tolerance for Runge-Kutta */\n", ode_epsabs);
		fprintf(fp, "ode_epsrel       \t%.17g\t/*relative error tolerance for Runge-Kutta */\n", ode_epsrel);
		fprintf(fp, "ode_nmax       \t%.2g\t/*maximum ode steps before reset */\n", (double)ode_nmax);
		fprintf(fp, "root_epsabs       \t%.17g\t/*absolute error tolerance for quasi-Newton */\n", root_epsabs);
		fprintf(fp, "root_epsrel       \t%.17g\t/*relative error tolerance for quasi-Newton */\n", root_epsrel);
		fprintf(fp, "alpha_tukey       \t%.8g\t/*value of parameter in Tukey window */\n", alpha_tukey);
	}
	else {
		printf("Failed to open file '%s'\n", parametersFilePathName);
	}
	if (fp != NULL) { fclose(fp); }		// COLM added to avoid errors
}


void readGlobalParameters(char *inFile) {
	paramFileBuffer = readParmetersFileToBuffer(inFile);
	VERBOSE = getIntParameterValueByName("Verbosity");
	getStringParameterValueByName("outputPath", SIM_DATA_OUTPUT);
	num_Threads = getIntParameterValueByName("numThreads");

	// Load in pulse parameters
	I_0 = getDoubleParameterValueByName("meanPumpIntensity");
	twoColorSH_amplitude = getDoubleParameterValueByName("twoColorRelativeIntensity");
	twoColorSH_phase = getDoubleParameterValueByName("twoColorPhase");
	tau = getDoubleParameterValueByName("pulseDuration");
	lambda_0 = getDoubleParameterValueByName("fundamentalWavelength"); 

	A_0 = sqrt(2.0 * I_0 / (epsilon_0*cLight));
	omega_0 = 2 * M_PI*cLight / lambda_0;

	// Load in domain parameters
	num_t = getIntParameterValueByName("numTimePoints");
	//num_t = pow(2, 17);
	domain_t = getDoubleParameterValueByName("timeDomainSize");

	//freqUpperCutoff = num_t / 2;
	freqLowerCutoff = getIntParameterValueByName("omegLowerCutoff");
	freqUpperCutoff = getIntParameterValueByName("omegUpperCutoff");
	numActiveOmega = num_t / 2 + 1;
	//numActiveOmega2 = numActiveOmega - (num_t / 2 + 1);
	l_0 = (num_t / 2 + 1)*(num_x / 2 + 1);

	sampleLayerThickness = getDoubleParameterValueByName("sampleLayerThickness");
	zStepMaterial1 = getDoubleParameterValueByName("initialZStep");
	LHSsourceLayerThickness = getDoubleParameterValueByName("LHSbufferThickness");
	RHSbufferLayerThickness = getDoubleParameterValueByName("RHSbufferThickness");
	
	alpha_tukey = getDoubleParameterValueByName("tukeyWindowAlpha");
}
