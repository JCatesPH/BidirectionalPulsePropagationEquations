// bragg_linear.cpp : Defines the entry point for the console application.
// ACMS version 1

//#include "stdafx.h"

#include "BPPE.h"
#include "Structure.h"
//#include "Utilities.h"
#include "createLayers.h"
#include <iomanip> 
#include <omp.h>
#include <vector>
#include <float.h>
#include <ctime>
#include <fenv.h>
#include <random>

//#define FFTW_WISDOM_TYPE FFTW_ESTIMATE
#define FFTW_WISDOM_TYPE FFTW_PATIENT
#define DO_DRUDE_MODEL
//#define DO_CONSTPLASMA

double dtime = omp_get_wtime();
using namespace std;

// Create GLOBAL structure and Material database
Structure myStructure("myFirstStructure");
Simulation mySimulation("andrewApplication1");
MaterialDB myMaterialsDB("myFirstMaterialDB");

// This block of vars were orignally inside main()
double zPosition = 0.0;
double sampleLayerThickness, I_0, A_0, tau, lambda_0, omega_0;
double twoColorSH_amplitude, twoColorSH_phase;
int num_t, freqLowerCutoff, freqUpperCutoff, numActiveOmega, numActiveOmega2, l_0, sizeRoot;
double domain_t, zStepMaterial1, alpha_tukey;
double *window;
vector<double> monitorZlocations;
int GSLerrorFlag, p, oFlag, VERBOSE;
double *omegaArray, *timeValuesArray, *kx, *ne, *y;
complex<double>*eFieldPlus, *eFieldMinus, *yp_init, *ym_init, *nl_k, *nl_p, *j_e;
fftw_plan nkForwardFFT, eFieldPlusForwardFFT, eFieldPlusBackwardFFT, eFieldMinusForwardFFT, eFieldMinusBackwardFFT, intBackwardFFT, npForwardFFT;
char *paramFileBuffer, SIM_DATA_OUTPUT[30];

int delmeFLAG = 0;

int main(int argc, char *argv[])
{	
	// Setting parameters using input file.
	paramFileBuffer = readParmetersFileToBuffer(argv[1]);
	VERBOSE = getIntParameterValueByName("Verbosity");
	getStringParameterValueByName("outputPath", SIM_DATA_OUTPUT);

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
	
	alpha_tukey = getDoubleParameterValueByName("tukeyWindowAlpha");
	//alpha_tukey = 0.0;

	pulseparam_type* sourceLeft = (pulseparam_type*)malloc(sizeof(pulseparam_type));
	pulseparam_type* sourceRight = (pulseparam_type*)malloc(sizeof(pulseparam_type));

	sourceLeft->A0 = sqrt(2.0 * I_0 / (epsilon_0*cLight));
	sourceLeft->omega0 = 2 * M_PI*cLight / lambda_0;
	sourceLeft->relativeIntensity = twoColorSH_amplitude;
	sourceLeft->relativePhase = twoColorSH_phase;
	sourceLeft->pulseDuration = tau;

	sourceRight->A0 = 0.0;
	sourceRight->omega0 = 2 * M_PI*cLight / lambda_0;
	sourceRight->relativeIntensity = twoColorSH_amplitude;
	sourceRight->relativePhase = twoColorSH_phase;
	sourceRight->pulseDuration = tau;
	//feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

	time_t now = time(0);
	char *datetime = ctime(&now);

	if (VERBOSE >= 0) {
		cout << "GBPPE code - ACMS Ver.0" << endl;
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


	generateLayers(myMaterialsDB, myStructure);
    setupPointMonitorLocations(myMaterialsDB, myStructure);
	/// VERY Early termination
	//printf("!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!: VERY Early Termination\n"); exit(-1);

	printf("Structure generated and point monitor locations set..\n");
	
	yp_init = (complex<double>*)malloc(sizeof(complex<double>) * numActiveOmega);
	ym_init = (complex<double>*)malloc(sizeof(complex<double>) * numActiveOmega);
	y = (double*)malloc(sizeof(double) * 4 * numActiveOmega);
	window = (double*)malloc(sizeof(double) * num_t);

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

	if (numDimensionsMinusOne == 1) {
		eFieldPlus = (complex<double>*)malloc(sizeof(complex<double>)*num_t*num_x);
		eFieldMinus = (complex<double>*)malloc(sizeof(complex<double>)*num_t*num_x);
		nl_k = (complex<double>*)malloc(sizeof(complex<double>)*num_t*num_x);
		nl_p = (complex<double>*)malloc(sizeof(complex<double>)*num_t*num_x);
		j_e = (complex<double>*)malloc(sizeof(complex<double>)*num_t*num_x);
		ne = (double*)malloc(sizeof(double)*num_t*num_x);
		kx = (double*)malloc(sizeof(double)*numActiveOmega);
		omegaArray = (double*)malloc(sizeof(double)*numActiveOmega);

		nkForwardFFT = fftw_plan_dft_2d(num_x, num_t, reinterpret_cast<fftw_complex*>(&nl_k[0]), reinterpret_cast<fftw_complex*>(&nl_k[0]), FFTW_FORWARD, FFTW_WISDOM_TYPE );
		eFieldPlusBackwardFFT = fftw_plan_dft_2d(num_x, num_t, reinterpret_cast<fftw_complex*>(&eFieldPlus[0]), reinterpret_cast<fftw_complex*>(&eFieldPlus[0]), FFTW_BACKWARD, FFTW_WISDOM_TYPE );
		eFieldPlusForwardFFT = fftw_plan_dft_2d(num_x, num_t, reinterpret_cast<fftw_complex*>(&eFieldPlus[0]), reinterpret_cast<fftw_complex*>(&eFieldPlus[0]), FFTW_FORWARD, FFTW_WISDOM_TYPE );
		npForwardFFT = fftw_plan_dft_2d(num_x, num_t, reinterpret_cast<fftw_complex*>(&j_e[0]), reinterpret_cast<fftw_complex*>(&nl_p[0]), FFTW_FORWARD, FFTW_WISDOM_TYPE );
		eFieldMinusBackwardFFT = fftw_plan_dft_2d(num_x, num_t, reinterpret_cast<fftw_complex*>(&eFieldMinus[0]), reinterpret_cast<fftw_complex*>(&eFieldMinus[0]), FFTW_BACKWARD, FFTW_WISDOM_TYPE );
		eFieldMinusForwardFFT = fftw_plan_dft_2d(num_x, num_t, reinterpret_cast<fftw_complex*>(&eFieldMinus[0]), reinterpret_cast<fftw_complex*>(&eFieldMinus[0]), FFTW_FORWARD, FFTW_WISDOM_TYPE );
		
	}
	else {
		eFieldPlus = (complex<double>*)malloc(sizeof(complex<double>)*num_t);
#ifdef WRITE_OUT_REFLECTANCE
		eFieldPlusBACKUPCOLM = (complex<double>*)malloc(sizeof(complex<double>) * num_t);
//		// DELETE THIS LOOP
//		double delMeSoon[12000];
//		double delMeSoon2[12000]; 
//#pragma loop(hint_parallel(4))
//		for (int i = 0; i < 12000; ++i)
//		{
//			delMeSoon[i] = delMeSoon2[i] * sqrt((double)i); // (y[i] + 1.0i * y[i]); // *exp(-1.0i * real(k[i]) * z);
//		}
#endif

		eFieldMinus = (complex<double>*)malloc(sizeof(complex<double>)*num_t);
		nl_k = (complex<double>*)malloc(sizeof(complex<double>)*num_t);
		nl_p = (complex<double>*)malloc(sizeof(complex<double>)*num_t);
		j_e = (complex<double>*)malloc(sizeof(complex<double>)*num_t);
		ne = (double*)malloc(sizeof(double)*num_t);
		omegaArray = (double*)malloc(sizeof(double)*num_t);
		kx = NULL;
	
		nkForwardFFT = fftw_plan_dft_1d(num_t, reinterpret_cast<fftw_complex*>(&nl_k[0]), reinterpret_cast<fftw_complex*>(&nl_k[0]), FFTW_FORWARD, FFTW_WISDOM_TYPE );
		eFieldPlusBackwardFFT = fftw_plan_dft_1d(num_t, reinterpret_cast<fftw_complex*>(&eFieldPlus[0]), reinterpret_cast<fftw_complex*>(&eFieldPlus[0]), FFTW_BACKWARD, FFTW_WISDOM_TYPE );
		eFieldPlusForwardFFT = fftw_plan_dft_1d(num_t, reinterpret_cast<fftw_complex*>(&eFieldPlus[0]), reinterpret_cast<fftw_complex*>(&eFieldPlus[0]), FFTW_FORWARD, FFTW_WISDOM_TYPE );
		npForwardFFT = fftw_plan_dft_1d(num_t, reinterpret_cast<fftw_complex*>(&j_e[0]), reinterpret_cast<fftw_complex*>(&nl_p[0]), FFTW_FORWARD, FFTW_WISDOM_TYPE );
		eFieldMinusBackwardFFT = fftw_plan_dft_1d(num_t, reinterpret_cast<fftw_complex*>(&eFieldMinus[0]), reinterpret_cast<fftw_complex*>(&eFieldMinus[0]), FFTW_BACKWARD, FFTW_WISDOM_TYPE );
		eFieldMinusForwardFFT = fftw_plan_dft_1d(num_t, reinterpret_cast<fftw_complex*>(&eFieldMinus[0]), reinterpret_cast<fftw_complex*>(&eFieldMinus[0]), FFTW_FORWARD, FFTW_WISDOM_TYPE );
		
	}


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


	fill_omg_k(omegaArray, kx);
	//DELME_ArgonDispersion(omegaArray);
	#ifdef DO_CONSTPLASMA
	DELME_AndrewPreformed(omegaArray);
	#endif
	createWindowFunc(alpha_tukey);
	//if (fp != NULL) { fclose(fp);  }
	
	if (VERBOSE >=3) printf("Generating the right-hand side source.\n");
	generateTwoColorPulse(eFieldMinus, eFieldMinusForwardFFT, ym_init, sourceRight);
	for (int i = 0; i < numActiveOmega; i++) {
		ym_init[i] = 0.0;
	} 
	if (VERBOSE >=3) printf("Generating the left-hand side source.\n");
    generateTwoColorPulse(eFieldPlus, eFieldPlusForwardFFT, yp_init, sourceLeft);
	writeInputSpectrum(yp_init);
	
	initializeY();
	doLinearProblem();

	std::string reldatpath = SIM_DATA_OUTPUT;
	myStructure.writeStructureLayoutToASCIIFile(reldatpath + "StructureLayout.txt");
	myStructure.writeStructureToDATFile(reldatpath + "Structure.dat");
	myMaterialsDB.writeMaterialDBToASCIIFile(reldatpath + "MaterialDatabase.txt");
	myStructure.writeBoundaryLayoutToASCIIFile(reldatpath + "BoundaryLayout.txt");

	/// Early termination
	//printf("############ WARNING ###########:  Early Termination\n"); exit(-1);
    printf("\n\n############ INFO ###########:  ENTERING THE NONLINEAR PART ##############################\n");

	iterateBPPE();

	dtime = omp_get_wtime() - dtime;
	if(VERBOSE >= 0) { cout << "Time in seconds is " << dtime << endl; }
	//cout << "Num threads set to  = " << omp_get_num_threads() << endl << endl;
	cout << endl << "Exiting program.." << endl << endl;

    return 0;
}

void setupPointMonitorLocations(MaterialDB& theMaterialDB, Structure& theStructure)
{

	//monitorZlocations.push_back(LHSsourceLayerThickness + 0.8e-6); 
	//monitorZlocations.push_back(LHSsourceLayerThickness + 2.0e-6); 
	//monitorZlocations.push_back(LHSsourceLayerThickness + 5.0e-6); 
	//monitorZlocations.push_back(LHSsourceLayerThickness + 10.0e-6); 
	//monitorZlocations.push_back(LHSsourceLayerThickness + 30.0e-6); 
	//monitorZlocations.push_back(LHSsourceLayerThickness + 10e-6); // 10 microns in plasma
	//monitorZlocations.push_back(LHSsourceLayerThickness + 20e-6); 
	//monitorZlocations.push_back(LHSsourceLayerThickness + 30e-6); 
	//monitorZlocations.push_back(LHSsourceLayerThickness + 40e-6); 
	//monitorZlocations.push_back(LHSsourceLayerThickness + 50e-6); 
	//monitorZlocations.push_back(LHSsourceLayerThickness + 75e-6);
	//monitorZlocations.push_back(LHSsourceLayerThickness + 100e-6); 

	//monitorZlocations.push_back(theStructure.getThickness() * 0.25);
	//monitorZlocations.push_back(theStructure.getThickness() * 0.5);
	//monitorZlocations.push_back(theStructure.getThickness() * 0.75);

	int num10ums = (myStructure.getThickness() - RHSbufferLayerThickness) / 10e-6;
	for (int n = 1; n < num10ums; n++){
		monitorZlocations.push_back(LHSsourceLayerThickness + n*10e-6);
	}
	monitorZlocations.push_back(myStructure.getThickness() - RHSbufferLayerThickness);

}

gsl_odeiv2_step * gslStep;
gsl_odeiv2_control * gslControl;
gsl_odeiv2_evolve * gslEvolve;

int mapU(const gsl_vector *ym_guess, void *rootparams, gsl_vector *f) {
    rootparam_type *rparams = reinterpret_cast<rootparam_type*>(rootparams);

	
	// Set the ODE params and set system up
    odeparam_type* params = fill_params(chi2_Material1, chi3_Material1, omegaArray, kx, 
			ne, j_e, myStructure.m_layers.begin()->getMaterial().getK(), eFieldPlus, eFieldMinus, nl_k, nl_p, 
			nkForwardFFT, eFieldPlusBackwardFFT, eFieldMinusBackwardFFT, npForwardFFT, myStructure.m_layers.front().getMaterial().getdoPlasmaCalc());
	
	//const gsl_odeiv2_step_type * stepType = gsl_odeiv2_step_rk4;
	//gsl_odeiv2_step * gslStep = gsl_odeiv2_step_alloc(stepType, 4 * numActiveOmega);
	//gsl_odeiv2_control * gslControl = gsl_odeiv2_control_y_new(ode_epsabs, ode_epsrel);
	//gsl_odeiv2_evolve * gslEvolve = gsl_odeiv2_evolve_alloc(4 * numActiveOmega);
	gsl_odeiv2_system sys = { func, NULL, (size_t)(4 * numActiveOmega), params };
	gsl_odeiv2_evolve_reset(gslEvolve);


    double nonlinear_time_initial, nonlinear_time;

	int numzReports, numZsteps = 0;
	double zRight, zStepSize;
	vector<double> endPoints;

    nonlinear_time_initial = omp_get_wtime();

	// Reset the left source and update guess
	for (int k = 0; k < numActiveOmega; k++){
		y[k] = real(yp_init[k]);
		y[k + numActiveOmega] = imag(yp_init[k]);
	}
	for (int k = 0; k < sizeRoot/2; k++){
		y[k + 2*numActiveOmega + freqLowerCutoff] = gsl_vector_get(ym_guess, k);
		y[k + 3*numActiveOmega + freqLowerCutoff] = gsl_vector_get(ym_guess, k + sizeRoot/2);
	}

	if (rparams->output == 1) {
		write_out_eFieldAndSpectrumAtZlocation(rparams->itnum, 0, y, 0.0, eFieldMinus, myStructure.m_layers.front().getMaterial().getK(), eFieldMinusBackwardFFT);
	}
    //cout << "  Going FORWARD through layers" << endl;
    for (std::list<Layer>::iterator lit = myStructure.m_layers.begin(); lit != myStructure.m_layers.end(); ++lit) {
        // Skip the LHS layer and the RHS layers
        if (lit->getLowSideBoundary() != NULL && lit->getHiSideBoundary() != NULL)
        {
            boundary(lit->getLowSideBoundary()->m_zPos, lit->getLowSideBoundary()->lowSideLayer->getMaterial().getK(), lit->getLowSideBoundary()->hiSideLayer->getMaterial().getK(), y);
            params->k = lit->getMaterial().getK();
            params->chi_2 = lit->getMaterial().getChi2();
            params->chi_3 = lit->getMaterial().getChi3();
            params->doPlasmaCalc = lit->getMaterial().getdoPlasmaCalc();
            params->mpi_k = lit->getMaterial().getmpi_k();
            params->mpi_sigmaK = lit->getMaterial().getmpi_sigmaK();
            params->ionE = lit->getMaterial().getIonizationEnergy();
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

                    GSLerrorFlag = gsl_odeiv2_evolve_apply(gslEvolve, gslControl, gslStep, &sys, &zPosition, zRight, &zStepSize, y);
                    
					if (GSLerrorFlag == GSL_SUCCESS) {
                        numZsteps++;
                    }
                    else {
                        printf("error: driver returned %d\n", GSLerrorFlag);
                        break;
                    }

                    nonlinear_time = omp_get_wtime() - nonlinear_time_initial;
                    if ((int)(nonlinear_time / 20) > numzReports) {
                        printf("  I = %d, step = %d, z = %.8g, t = %d s\n", rparams->itnum, numZsteps, zPosition, (int)nonlinear_time);
                        numzReports++;
                    }
                }

				if (rparams->output == 1) {
					//printf("  Outputting Point Monitor file at Z location %d[nm]... \n", (int)round(zPosition * 1.0e9));
					write_multicolumnMonitor(rparams->itnum, zPosition, y, params);

				}
                
            }

            endPoints.clear();

            //integrate(zPosition, zStepSize, params, y, integral);
            
        }
        //  Finally do the lowside boundary of the last assumed Vacuum layer
        if (lit->getHiSideBoundary() == NULL) {
            boundary(lit->getLowSideBoundary()->m_zPos, lit->getLowSideBoundary()->lowSideLayer->getMaterial().getK(), lit->getLowSideBoundary()->hiSideLayer->getMaterial().getK(), y);
        }
    }

	for (int k = 0; k < sizeRoot/2; k++){
		gsl_vector_set(f, k, y[k + 2 * numActiveOmega + freqLowerCutoff] - real(ym_init[k + freqLowerCutoff]));
		gsl_vector_set(f, k + sizeRoot/2, y[k + 3 * numActiveOmega + freqLowerCutoff] - imag(ym_init[k + freqLowerCutoff]));
	}


	if (rparams->output == 1) {
		write_out_eFieldAndSpectrumAtZlocation(rparams->itnum, 1, y, myStructure.getThickness(), eFieldPlus, myStructure.m_layers.back().getMaterial().getK(), eFieldPlusBackwardFFT);
	}
    
	nonlinear_time = omp_get_wtime() - nonlinear_time_initial;
	//if (rparams->output == 1) {
		//printf("Iteration %d completed in %.2f seconds with %d steps.\n", rparams->itnum, nonlinear_time, numZsteps);
		//cout << "Iteration " << rparams->itnum <<  " completed in " <<  nonlinear_time << "seconds with" << numZsteps << "steps." << endl;
		//fflush(stdout);
	//}

	for (int k = 0; k < numActiveOmega; k++){
		y[k + 2*numActiveOmega] = real(ym_init[k]);
		y[k + 3*numActiveOmega] = imag(ym_init[k]);
	}

	myStructure.doBackwardPassThroughAllBoundaries(y);


	//gsl_odeiv2_control_free(gslControl);
    //gsl_odeiv2_evolve_free(gslEvolve);
    //gsl_odeiv2_step_free(gslStep);
    
    return GSL_SUCCESS;
}

void iterateBPPE()
{
	// Find the size of the problem with omega cutoffs
	sizeRoot = 2*(freqUpperCutoff - freqLowerCutoff + 1);

	// Initialize param struct for root solver
    rootparam_type *rparams = (rootparam_type*)malloc(sizeof(rootparam_type));
    rparams->itnum = 1;
	rparams->output = 0;

	// Initialize time and status variables
	int status;
	double nonlinear_time_initial, nonlinear_time, nonlinear_time_tmp, nonlinear_time_total;
	nonlinear_time_initial = omp_get_wtime();
	nonlinear_time_tmp = nonlinear_time_initial;

	// Initialize multiroot objects
	printf("Allocating multiroot solver\n");
	const gsl_multiroot_fsolver_type *T;
    gsl_multiroot_fsolver *s;
	T = gsl_multiroot_fsolver_hybrids;
	s = gsl_multiroot_fsolver_alloc(T, sizeRoot);
	
	printf("Allocating initial guess\n");
	gsl_vector *u = gsl_vector_alloc(sizeRoot);
	//u->data = y;

	// Create pseudo-random number generator for initial guess
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<double> dis(-1.0, 1.0);
	// Set the initial guess with y and uniform r.v.
	for (int k = 0; k < sizeRoot/2; k++){
		gsl_vector_set(u, k, y[k + 2*numActiveOmega + freqLowerCutoff] + dis(gen) * INITIAL_GUESS_SEED_VALUE);
		gsl_vector_set(u, k + sizeRoot/2, y[k + 3*numActiveOmega + freqLowerCutoff] + dis(gen) * INITIAL_GUESS_SEED_VALUE);
	}

	// Initialize GSL ODE objects
	const gsl_odeiv2_step_type * stepType = gsl_odeiv2_step_rk4;
	gslStep = gsl_odeiv2_step_alloc(stepType, 4 * numActiveOmega);
	gslControl = gsl_odeiv2_control_y_new(ode_epsabs, ode_epsrel);
	gslEvolve = gsl_odeiv2_evolve_alloc(4 * numActiveOmega);

	// Tell GSL multiroot the function and initial guess
	printf("Setting multiroot function\n");
	gsl_multiroot_function f = {&mapU, sizeRoot, rparams};
	gsl_multiroot_fsolver_set(s, &f, u);
	printf("Finished setting multiroot function\n");

	rparams->output = 1;
	do
	{
		printf("Starting iteration %d\n", rparams->itnum);
		fflush(stdout);
		nonlinear_time_tmp = omp_get_wtime();
		status = gsl_multiroot_fsolver_iterate(s);

		if (status) break;

		//status = gsl_multiroot_test_residual(s->f, root_epsabs);
		status = gsl_multiroot_test_delta(s->dx, s->x, root_epsabs, root_epsrel);

		nonlinear_time = omp_get_wtime() - nonlinear_time_tmp;
		printf("Iteration %d completed in %.2f seconds.\n", rparams->itnum, nonlinear_time);
		rparams->itnum = rparams->itnum + 1;
		fflush(stdout);
	}
	while (status == GSL_CONTINUE && rparams->itnum < 25000);

	
	nonlinear_time_total = omp_get_wtime() - nonlinear_time_initial;
	printf("  Multiroot solver completed in %.2f seconds.\n\n", nonlinear_time_total);

	// Run the map one last time to output spectra
	printf("Performing final iteration with output enabled..\n");
	rparams->output = 1;
	gsl_vector *tmpf = gsl_vector_alloc(2*numActiveOmega);
	mapU(s->x, rparams, tmpf);

	// Freeing ODE memory
	gsl_odeiv2_control_free(gslControl);
    gsl_odeiv2_evolve_free(gslEvolve);
    gsl_odeiv2_step_free(gslStep);

	// Free solver memory
	printf("Freeing solver memory.\n");
    gsl_multiroot_fsolver_free(s);
	gsl_vector_free(u);
	gsl_vector_free(tmpf);


}


void fill_omg_k(double*omg, double*kx) {

		for (int i = 0; i < num_t; i++)
		{
			if (i <= num_t / 2) {
				omg[i] = (M_PI / domain_t)*i;
			}
			else {
				omg[i] = (M_PI / domain_t)*((double)i - num_t);
			}
		}
		// COLM NEW this line replaces following comment block
		myMaterialsDB.initAllMaterialKs(omg, numActiveOmega);

	return;
}

#define ANDREW_PREFORMED_METHOD1
void DELME_AndrewPreformed(double* omg) {
	Material *plasmaMat;
	plasmaMat = myMaterialsDB.getMaterialByName("PlasmaMat");

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

	plasmaMat->m_k[0] = 0.0;
	fprintf(fp, "%.7g\t%.17g\t%.17g \n", 0.0, real(plasmaMat->m_k[0]), imag(plasmaMat->m_k[0]));
	for (int i = 1; i < numActiveOmega; i++)
	{
		complex<double> k0 = omg[i] / cLight;
		//plasmaMat->m_k[i] = k0 * sqrt(complex<double>(1.0 - pow(omega_plasma / omg[i], 2)));
		plasmaMat->m_k[i] = k0 * sqrt(1.0 - pow(omega_plasma, 2) / (pow(omg[i], 2) + 1.0i * omg[i] / tauCollision));
		fprintf(fp, "%.7g\t%.17g\t%.17g \n", omg[i], real(plasmaMat->m_k[i]), imag(plasmaMat->m_k[i]));
	}
	if (fp != NULL) { fclose(fp); }
}

void DELME_ArgonDispersion(double* omg) {
	Material *ArgonMat;
	complex<double> n0;
	double lambda;
	double raylighFactor = 32.0 * pow(M_PI, 3) / (3.0 * pow(num_atoms,2)) * num_atoms;
	double rayleighScattering = 0.0;
	ArgonMat = myMaterialsDB.getMaterialByName("Argon");

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
			n0 = 1.0 + 2.50141e-3/(91.012-pow(lambda,-2)) + 5.00283e-4/(87.892-pow(lambda,-2)) + 5.22343e-2/(214.02-pow(lambda,-2));
		}
		else {
			n0 = 1.0 + 6.7867e-5 + 3.0182943e-2/(144.0-pow(lambda,-2)); // From Peck and Fisher (1964), valid in range: 0.47-2.06 um (0 deg C)
		}

		// Fitted dispersion (in Mathematica)
		n0 = 1.00026 + 1.37647e-5 / lambda - 1.86046e-6 / pow(lambda,2) + 3.78465e-7 / pow(lambda,3);
		
		//n0 += 1.0i * 8.0 * pow(M_PI, 3) / (3 * num_atoms * pow(lambda,4)) * pow(pow(n0, 2) - 1.0, 2);
		rayleighScattering = raylighFactor / pow(lambda*1e-6, 4) * pow(real(n0)-1.0, 2);
		//n0 += 1.0i * rayleighScattering;

		ArgonMat->m_k[i] = omg[i] * n0 / cLight;

		fprintf(fp, "%.7g\t%.17g\t%.17g \n", lambda, real(n0), imag(n0));
	}
	if (fp != NULL) { fclose(fp); }
}

//void copy_omg_k_ToMaterial(MaterialDB aMaterialDB, double* kx, complex<double>* k_0, complex<double>* k_1, complex<double>* k_2, complex<double>* k_3) {
//
//}
void initializeY()
{
	for (int i = 0; i < numActiveOmega; i++)
	{
		y[i] = real(yp_init[i]);
		y[i + numActiveOmega] = imag(yp_init[i]);
		y[i + 2*numActiveOmega] = 0.0;
		y[i + 3*numActiveOmega] = 0.0;
	}

	for (int i = 0; i < freqLowerCutoff; i++) {
		y[i] = 0.0;
		y[i + numActiveOmega] = 0.0;
		y[i + 2*numActiveOmega] = 0.0;
		y[i + 3*numActiveOmega] = 0.0;
	}
	for (int i = freqUpperCutoff; i < numActiveOmega; i++) {
		y[i] = 0.0;
		y[i + numActiveOmega] = 0.0;
		y[i + 2*numActiveOmega] = 0.0;
		y[i + 3*numActiveOmega] = 0.0;
	}

}

void doLinearProblem() {
	for (int its = 0; its < 3; its++) {
		// Pass forward through boundaries 
		myStructure.doForwardPassThroughAllBoundaries(y);

		// Set right source
		for (int i = 0; i < numActiveOmega; i++)
		{
			y[i + 2*numActiveOmega] = real(ym_init[i]);
			y[i + 3*numActiveOmega] = imag(ym_init[i]);
		}

		// Pass backward through boundaries
		myStructure.doBackwardPassThroughAllBoundaries(y);

		// Reset left source to ensure consistency
		for (int i = 0; i < numActiveOmega; i++)
		{
			y[i] = real(yp_init[i]);
			y[i + numActiveOmega] = imag(yp_init[i]);
		}
	}
	
	// Output the reflected spectrum
	write_out_eFieldAndSpectrumAtZlocation(0, 0, y, 0.0, eFieldMinus, myStructure.m_layers.front().getMaterial().getK(), eFieldMinusBackwardFFT);

	// Pass forward through boundaries 
	myStructure.doForwardPassThroughAllBoundaries(y);

	// Output the transmitted spectrum 
	write_out_eFieldAndSpectrumAtZlocation(0, 1, y, myStructure.getThickness(), eFieldPlus, myStructure.m_layers.back().getMaterial().getK(), eFieldPlusBackwardFFT);

	// Set right source
	/* for (int i = 0; i < numActiveOmega; i++)
	{
		y[i + 2*numActiveOmega] = real(ym_init[i]);
		y[i + 3*numActiveOmega] = imag(ym_init[i]);
	} */

	// Pass backward through boundaries
	myStructure.doBackwardPassThroughAllBoundaries(y);

	// Reset left source to ensure consistency
	for (int i = 0; i < numActiveOmega; i++)
	{
		y[i] = real(yp_init[i]);
		y[i + numActiveOmega] = imag(yp_init[i]);
	}

}

void writeInputEfield(std::complex<double>* ee_p)
{
	char inputEfieldFilePathName[STRING_BUFFER_SIZE];
	snprintf(inputEfieldFilePathName, sizeof(char) * STRING_BUFFER_SIZE, "%sInputEfield_1D.dat", SIM_DATA_OUTPUT);
	FILE* fp;
	fp = fopen(inputEfieldFilePathName, "w");
	double dt = domain_t / double(num_t);

	if (fp != NULL)
	{
		fprintf(fp, "# time [sec]\tEfield [V/m]\n");
		for (int i = 0; i < num_t; i++)
		{
			//fprintf(fp, "%.10lf \n", real(eFieldPlus[i]));							//PARIS	formated in single column file
			fprintf(fp, "%.7g\t%.17g \n", i * dt, real(ee_p[i]));   // COLM export in one column format
		}
	}
	else {
		printf("Failed to open file '%s'\n", inputEfieldFilePathName);

		exit(-1);
	}
	if (fp != NULL) { fclose(fp); }
}

void writeInputSpectrum(std::complex<double>* y_init)
{
	char inputEfieldFilePathName2[STRING_BUFFER_SIZE];
	snprintf(inputEfieldFilePathName2, sizeof(char) * STRING_BUFFER_SIZE, "%sInputSpectrum_1D.dat", SIM_DATA_OUTPUT);

	//printf("Input Spectrum file: %s\n", inputEfieldFilePathName2);

	FILE* fp_spectrum;
	fp_spectrum = fopen(inputEfieldFilePathName2, "w");
	if (fp_spectrum != NULL)
	{
		fprintf(fp_spectrum, "# omega [Hz]\tReal []\tImag []\tAbs []\n");
		for (int i = 0; i < numActiveOmega; i++)
		//for (int i = 0; i < num_t; i++)
		{
			//fprintf(fp, "%.10lf \n", real(eFieldPlus[i]));							//PARIS	formated in single column file
			//fprintf(fp, "%.7g\t%.17g \n", i * domain_t / num_t, real(ee_p[i]));		// COLM export in one column format
			//fprintf(fp_spectrum, "%g \t %+.17g \t %+.17g\n", (M_PI / domain_t)* i, real(yp_init[i]), imag(yp_init[i]));		// COLM input spectrum output
			fprintf(fp_spectrum, "%g \t %+.17g \t %+.17g \t %+.17g\n", (2 * M_PI / domain_t) * i, real(y_init[i]), imag(y_init[i]), abs(y_init[i]));
			// COLM outputing abs(Sp) too
#ifdef WRITE_OUT_REFLECTANCE																																						// COLM make a backup of input spectrum to calculate reflectance later on
			eFieldPlusBACKUPCOLM[i] = yp_init[i];
#endif
		}
	}
	else {
		printf("Failed to open file '%s'\n", inputEfieldFilePathName2);

		exit(-1);
	}
	if (fp_spectrum != NULL) { fclose(fp_spectrum); }
}

void createWindowFunc(double alpha){
	// Implements a Tukey window
	for (int i = 0; i < num_t/8; i++) {
		window[i] = 0.0;
	}
	for (int i=num_t/8; i < (int)((1+alpha)*num_t/8); i++) {
		window[i] = 0.5 * (1.0 - cos(4*M_PI*(i-num_t/4)/(alpha*num_t)));
	}
	for (int i=(int)((1+alpha)*num_t/8); i <= num_t/2; i++) {
		window[i] = 1.0;
	}
	for (int i=num_t/8; i <= num_t/2; i++) {
		window[num_t-i] = window[i];
	}
	for (int i = 7*num_t/8; i < num_t; i++) {
		window[i] = 0.0;
	}
	/* for (int i=0; i < (int)(alpha*num_t/2); i++) {
		window[i] = 0.5 * (1.0 - cos(2.0*M_PI*i/(alpha*num_t))) + 1e-15;
	}
	for (int i=(int)(alpha*num_t/2); i <= num_t/2; i++) {
		window[i] = 1.0 + 1e-15;
	}
	for (int i=1; i <= num_t/2; i++) {
		window[num_t-i] = window[i];
	} */

	char windowFile[STRING_BUFFER_SIZE];
	snprintf(windowFile, sizeof(char) * STRING_BUFFER_SIZE, "%swindowFunc.dat", SIM_DATA_OUTPUT);
	FILE* fp;
	fp = fopen(windowFile, "w");
	for (int i=0; i < num_t; i++){
		fprintf(fp, "%.16f \t %.16f\n", i*domain_t/num_t, window[i]);
	}
	fclose(fp);
}

void createWindowFunc(){
	// Implements a Hann/Hamming window
	double a0 = 25.0/46.0;
	for (int i = 0; i < num_t; i++) {
		window[i] = a0 - (1.0 - a0) * cos(2.0*M_PI*i/num_t);
	}

	char windowFile[STRING_BUFFER_SIZE];
	snprintf(windowFile, sizeof(char) * STRING_BUFFER_SIZE, "%swindowFunc.dat", SIM_DATA_OUTPUT);
	FILE* fp;
	fp = fopen(windowFile, "w");
	for (int i=0; i < num_t; i++){
		fprintf(fp, "%.16f \t %.16f\n", i*domain_t/num_t, window[i]);
	}
	fclose(fp);
}


void normalizeFFT(complex<double>* arr, int type) {
	
	double c = sqrt((double)num_t);
	if (type == 1) {
		for (int i=0; i < num_t; i++) {
			arr[i] = arr[i] / c;
		}
	}
	else if (type == 2) {
		for (int i=0; i < num_t; i++) {
			arr[i] = arr[i] * window[i] / c;
		}
	}
	else {
		printf("WARNING: Invalid option given to FFT normalization function.\n");
		c = (double(num_t));
		for (int i=0; i < num_t; i++) {
			arr[i] = arr[i] / c;
		}
	}

}

void applyWindow(complex<double>* arr){
	for (int i=0; i < num_t; i++) {
		arr[i] = arr[i] * window[i];
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
    normalizeFFT(ee, 1);

	for (int i = 0; i < numActiveOmega; i++) {
		source[i] = ee[i];
	}

}



odeparam_type* fill_params(double chi_2, double chi_3, double* omg, double* kx, double* ne, complex<double>* j_e, complex<double>* k, complex<double>* ee_p, complex<double>* ee_m, complex<double>* nl_k, complex<double>* nl_p, fftw_plan nk_f, fftw_plan ep_b, fftw_plan em_b, fftw_plan np_f, int plasmaBool) {

	odeparam_type* r = (odeparam_type*)malloc(sizeof(odeparam_type));
	if (r != NULL)
	{
		r->chi_2 = chi_2;
		r->chi_3 = chi_3;
		r->omega = omg;
		r->kx = kx;
		r->rho = ne;
		r->j_e = j_e;
		r->k = k;
		r->ee_p = ee_p;
		r->ee_m = ee_m;
		r->nl_k = nl_k;
		r->nl_p = nl_p;
		r->nk_f = nk_f;
		r->ep_b = ep_b;
		r->em_b = em_b;
		r->np_f = np_f;
	}
	else {
		printf("Could not maloc space in fill_params()\n");
		exit(-1);
	}
	return r;
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
		fprintf(fp, "cLight         \t%.17g\t/*speed of light */\n", cLight);
		fprintf(fp, "epsilon_0      \t%.17g\t/*permittivity of free space */\n", epsilon_0);
		fprintf(fp, "Znaught        \t%.17g\t/*impedance of free space */\n", Znaught);
		fprintf(fp, "charge_e       \t%.17g\t/*electron charge */\n", charge_e);
		fprintf(fp, "mass_e         \t%.17g\t/*electron mass */\n", mass_e);
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
		fprintf(fp, "waist_x         \t%.17g\t/*pulse-waist fwhm */\n", waist_x);
		fprintf(fp, "num_t           \t%.17g\t/*number of time pos */\n", (double)num_t);
		fprintf(fp, "num_x           \t%.17g\t/*number of x pos */\n", (double)num_x);
		fprintf(fp, "domain_t        \t%.17g\t/*time domain */\n", domain_t);
		fprintf(fp, "domain_x        \t%.17g\t/*x domain */\n", domain_x);
		fprintf(fp, "num_atoms       \t%.17g\t/*number of atoms in gas */\n", num_atoms);
		fprintf(fp, "rho_0           \t%.17g\t/* initial electron density */\n", rho_0);
		fprintf(fp, "j_e0            \t%.17g\t/* initial current density */\n", j_e0);
		fprintf(fp, "freqUpperCutoff \t%d\t/* upper frequency cut-off */\n", freqUpperCutoff);
		fprintf(fp, "freqLowerCutoff \t%d\t/*lower frequency cut-off */\n", freqLowerCutoff);
		fprintf(fp, "num_iterations  \t%d\t/*number of BPPE iterations*/ */\n", num_iterations);
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
		fprintf(fp, "epsabs       \t%.17g\t/*absolute error tolerance for Runge-Kutta */\n", ode_epsabs);
		fprintf(fp, "epsrel       \t%.17g\t/*relative error tolerance for Runge-Kutta */\n", ode_epsrel);
		fprintf(fp, "ode_nmax       \t%.2g\t/*maximum ode steps before reset */\n", (double)ode_nmax);

		fprintf(fp, "alpha_tukey       \t%.8g\t/*value of parameter in Tukey window */\n", alpha_tukey);
	}
	else {
		printf("Failed to open file '%s'\n", parametersFilePathName);
	}
	if (fp != NULL) { fclose(fp); }		// COLM added to avoid errors
}



void write_out_eFieldAndSpectrumAtZlocation(int num, int j, double*y, double z, complex<double>*ee, complex<double>*k, fftw_plan e_b) {

	char efieldFilePathName[STRING_BUFFER_SIZE];
	char spectrumFilePathName[STRING_BUFFER_SIZE];
	// COLMTODO E_ij -> E_iter_i_Reflect(0) ,Transmitted(1)
	snprintf(efieldFilePathName, sizeof(char) * STRING_BUFFER_SIZE, "%sEfield_iteration_%i_%s.dat", SIM_DATA_OUTPUT, num, (j == 0 ? "Reflected" : "Transmitted"));
	snprintf(spectrumFilePathName, sizeof(char) * STRING_BUFFER_SIZE, "%sSpectrum_iteration_%i_%s.dat", SIM_DATA_OUTPUT, num, (j == 0 ? "Reflected" : "Transmitted"));
	//snprintf(buffer2, sizeof(char) * STRING_BUFFER_SIZE, "Snew_%i%i.dat", Iteration_number, j);

	if 	(z>myStructure.getThickness()) {
		printf("Request to write field at point outside structure %g  Structure range 0<->%g\n", z, myStructure.getThickness());
	}


	if (j == 1) {

		for (int i = 0; i <= num_t / 2; i++)
		{
			//ee[i] = (y[i] + 1.0i * y[i + num_t / 2 + 1]) * exp(-1.0i * k[i] * (zPosition + RHSbufferLayerThickness));		// phase corrections by Andrew (2021-01-25)
			//ee[i] = (y[i] + 1.0i * y[i + num_t / 2 + 1]) * exp(1.0i * real(k[i]) * z) * exp(-1.0 * abs(imag(k[i])) * z);
			ee[i] = (y[i] + 1.0i * y[i + num_t / 2 + 1]) * exp(1.0i * k[i] * z);
		}

		for (int i = 1; i < num_t / 2; i++)
		{
			//ee[num_t - i] = (y[i] - 1.0i*y[i + num_t / 2 + 1])*exp(1.0i*conj(k[i])*(zPosition + RHSbufferLayerThickness));		// phase corrections by Andrew (2021-01-25)
			//ee[num_t - i] = (y[i] - 1.0i * y[i + num_t / 2 + 1]) * exp(-1.0i * real(k[i]) * z) * exp(-1.0 * abs(imag(k[i])) * z);
			ee[num_t - i] = (y[i] - 1.0i * y[i + num_t / 2 + 1]) * exp(-1.0i * k[i] * z);
		}
	}
	else {
		for (int i = 0; i <= num_t / 2; i++)
		{
			ee[i] = (y[i + num_t + 2] + 1.0i*y[i + 3 * num_t / 2 + 3]);
		}

		for (int i = 1; i < num_t / 2; i++)
		{
			ee[num_t - i] = (y[i + num_t + 2] - 1.0i*y[i + 3 * num_t / 2 + 3]);
		}
	}

	// PARIS output spectrum
	FILE* fp2;
	//errno_t err2;
	fp2 = fopen(spectrumFilePathName, "w");
	if (fp2 != NULL)
	{
		fprintf(fp2, "# omega [Hz]\tReal []\tImag []\tAbs []      recorded at zPosition=%g[microns]\n", z * 1e6);

		for (int i = 0; i < num_t; i++)
		{
			//fprintf(fp2, "%g \t %+.17g \t %+.17g\t %+.17g\n", (M_PI/domain_t)*i, real(ee[i]), imag(ee[i]), abs(ee[i]));		// COLM export in one column format
			if (i <= num_t / 2) {
				fprintf(fp2, "%g \t %+.17g \t %+.17g \t %+.17g\n", (2 * M_PI / domain_t) * i, real(ee[i]), imag(ee[i]), abs(ee[i]));
			}
			else {
				fprintf(fp2, "%g \t %+.17g \t %+.17g \t %+.17g\n", (2 * M_PI / domain_t) * (i - num_t), real(ee[i]), imag(ee[i]), abs(ee[i]));
			}
		}
	}
	else {
		printf("Failed to open file '%s'\n", spectrumFilePathName);
	}
	if(fp2 != NULL) { fclose(fp2);  }		// COLM added to avoid errors

#ifdef WRITE_OUT_REFLECTANCE
	if (j == 0)
	{
		// COLM output Reflectance spectrum
		char reflectanceFilePathName[STRING_BUFFER_SIZE];
		snprintf(reflectanceFilePathName, sizeof(char) * STRING_BUFFER_SIZE, "%sReflectanceSpectrum_iteration_%i_%s.dat", SIM_DATA_OUTPUT, num, (j == 0 ? "Reflected" : "Transmitted"));
		FILE* fp3;
		//errno_t err3;
		err3 = fopen_s(&fp3, reflectanceFilePathName, "w");
		if (fp3 != NULL)
		{
			fprintf(fp3, "# omega [Hz] \tReflectance []\n");

			for (int i = 0; i <= numActiveOmega; i++)
			{
				fprintf(fp3, "%g \t %+.17g \n", (M_PI / domain_t) * i, abs((ee[i] * ee[i]) / (eFieldPlusBACKUPCOLM[i] * eFieldPlusBACKUPCOLM[i])));		// COLM export in one column format
			}
		}
		else {
			printf("Failed to open file '%s'\n", reflectanceFilePathName);
		}
		if (fp3 != NULL) { fclose(fp3); }		// COLM added to avoid errors 
	}
#endif

	fftw_execute(e_b);
	normalizeFFT(ee, 2);

	FILE *fp;
	//errno_t err;
	fp = fopen(efieldFilePathName, "w");
	if (fp != NULL)
	{
		fprintf(fp, "# time [sec]\tEfield [V/m]    recorded at zPosition=%g[microns]\n",z*1e6);
		for (int i = 0; i < num_t; i++)
		{
			fprintf(fp, "%.7g\t%.17g \n", i*domain_t/num_t, real(ee[i]));		// COLM export in one column format
		}
	}
	else {
		printf("Failed to open file '%s'\n", efieldFilePathName);
	}
	if (fp != NULL) { fclose(fp);  }		// COLM added to avoid errors
	return;
}



int func(double z, const double y[], double f[], void *odep) {

	odeparam_type *p = reinterpret_cast<odeparam_type*>(odep);

	const int num_tOver2 = num_t / 2;
	const double clightSquared = pow(cLight, 2);
	const double num_td = (double)num_t;

	#pragma omp parallel for
	for (int i = 0; i <= num_tOver2; i++)
	{
		//const complex<double> phaseFactor = exp(1.0i * real(p->k[i]) * z) * exp(-1.0 * abs(imag(p->k[i])) * z);
		const complex<double> phaseFactor = exp(1.0i * p->k[i] * z);
		p->ee_p[i] = (y[i] + 1.0i * y[i + num_tOver2 + 1]) * phaseFactor;
		p->ee_m[num_t - i] = (y[i + num_t + 2] - 1.0i * y[i + 3 * num_tOver2 + 3]) * phaseFactor;

		if (i > 0 && i < num_tOver2) {
			//const complex<double> phaseFactor2 = exp(-1.0i * real(p->k[i]) * z) * exp(-1.0 * abs(imag(p->k[i])) * z);
			const complex<double> phaseFactor2 = exp(-1.0i * p->k[i] * z);
			p->ee_p[num_t - i] = (y[i] - 1.0i * y[i + num_tOver2 + 1]) * phaseFactor2;
			p->ee_m[i] = (y[i + num_t + 2] + 1.0i * y[i + 3 * num_tOver2 + 3]) * phaseFactor2;
		}
		
	}

	
	// DELETE ME : Drude model
	#ifdef DO_DRUDE_MODEL
	const double sig0 = rho_0 * pow(charge_e, 2) * tauCollision / mass_e;
	const double omeg_p2 = rho_0 * pow(charge_e, 2) / (mass_e * epsilon_0); 
	#pragma omp parallel for
	for (int i = 0; i < num_t; i++)
	{	
		// Calculate the current density
		//p->nl_p[i] = sig0 / (1.0 - 1.0i * p->omega[i] * tauCollision) * (p->ee_p[i]+p->ee_m[i]);
		//const double w = sig0 / (1.0 + pow(p->omega[i]*tauCollision, 2));
		//p->nl_p[i] = w * (1.0 + 1.0i * p->omega[i] * tauCollision) * (p->ee_p[i]+p->ee_m[i]);
		p->nl_p[i] = 0.0;

		// Calculate the polarization
		//p->nl_k[i] = (1.0 - omeg_p2 / (pow(p->omega[i], 2) + 1.0i * p->omega[i] / tauCollision)) * (p->ee_p[i]+p->ee_m[i]);
		p->nl_k[i] = (-omeg_p2 / (pow(p->omega[i], 2) + 1.0i * p->omega[i] / tauCollision)) * (p->ee_p[i] + p->ee_m[i]);
		//p->nl_k[i] = p->nl_k[i] / sqrt((double)num_t);
		//p->nl_k[i] = p->nl_k[i] / sqrt(2.0*M_PI);
	}

	//fftw_execute(p->nk_f);
	//normalizeFFT(nl_k, 1);

	// ---------------------------------
	#else

	fftw_execute(p->ep_b);
	fftw_execute(p->em_b);
	normalizeFFT(p->ee_p, 1);
	normalizeFFT(p->ee_m, 1);

	#pragma omp parallel for
	for (int i = 0; i < num_t; i++)
	{
		//p->ee_p[i] = p->ee_p[i] / num_td;
		//p->ee_m[i] = p->ee_m[i] / num_td;
		p->nl_k[i] = p->chi_2 * pow(real(p->ee_p[i] + p->ee_m[i]), 2) + p->chi_3 * pow(real(p->ee_p[i] + p->ee_m[i]), 3);
	}
	
	fftw_execute(p->nk_f);
	normalizeFFT(nl_k, 1);

	
	if (p->doPlasmaCalc == 2) {

		// POSSIBLE ERROR WHY FACTOR 2.0 in following ht calculation???
		double ht = domain_t / num_td;
		double neutrals = num_atoms - rho_0;                          // Neutral particles
		double electrons = rho_0;                      // background Electrons
		double change = 0.0;    
		double current = 0.0;                                       // Current to be exported to UPPE
		double current_change = 0.0e0;                              // Current change for differential equation
		double ve = 0.0; //1/tauCollision;
		double fv1 = 0.0e0, fv2 = 0.0e0;

		p->rho[0] = electrons;
		for (int i = 0; i < num_t - 1; i++)
		{
			double fieldIntensity = ( pow(real(p->ee_p[i] + p->ee_m[i]),2) ) / Znaught ;  // MIRO real+real
			change = neutrals * p->mpi_sigmaK * ht * pow(fieldIntensity, p->mpi_k);
			electrons += change;
			neutrals -= change;
			// Don't allow neutrals dip below zero
			if (neutrals < 0.0) neutrals = 0.0;
			if (electrons > num_atoms) electrons = num_atoms;

			p->rho[i + 1] = electrons;
		}

		p->j_e[0] = j_e0;
		fv1 = j_e0;
		for (int i = 0; i < num_t - 1; i++)
		{
			//p->j_e[i + 1] = (1.0 - ht / tauCollision)*p->j_e[i] + ht * pow(charge_e, 2) / mass_e * p->rho[i] * real(p->eFieldPlus[i] + p->eFieldMinus[i]);
			// Andrew original
			//p->j_e[i + 1] = p->j_e[i] * exp(-1.0 * ht / tauCollision) + pow(charge_e, 2) / mass_e * 0.5 * ht * (p->rho[i] * real(p->ee_p[i] + p->ee_m[i]) * exp(-1.0 * ht / tauCollision) + p->rho[i + 1] * real(p->ee_p[i + 1] + p->ee_m[i + 1]));
			
			
			// CALCULATING THE CURRENT for UPPE (based on Ewan's notes Jan-24 2018)
			//  First order method	
			//	current_change = delta_t*(pow(cnst_e,2)/cnst_me)*electrons*real(aptr[t])-delta_t*ve*current;	
			// Second order method (Kolja)
			fv2 = real(p->ee_p[i + 1] + p->ee_m[i + 1]);
			current_change = ht * (pow(charge_e, 2) / mass_e) * electrons * (fv1 + fv2) * 0.5 - ht * ve * current;
			fv1 = fv2;
			current += current_change;
			p->j_e[i + 1] = current;
		}

		//if (z == distanceSourceToSample)
		

		// I THINK.. this takes the Current j_e and forewardFFTs it into the array called nl_p which is then used in the Integrate()
		//applyWindow(p->j_e);
		fftw_execute(p->np_f);
		normalizeFFT(nl_p, 1);

	}
	else if (p->doPlasmaCalc == 1){
		double ht = domain_t / num_td;
		double neutrals = num_atoms - rho_0;  // Neutral particles
		double electrons = rho_0; 
		double change = 0.0;    
		double current = 0.0;                                       // Current to be exported to UPPE
		double current_change = 0.0e0;                              // Current change for differential equation
		double ve = 1/tauCollision;
		double fv1 = 0.0e0, fv2 = 0.0e0;
		const double nu_a = 4.13e16; // [Hz]
		const double E_a = 5.14e11; // [V/m]
		const double U_H = 13.6; // Ionization potential of hydrogen [eV]
		double potentialFrac = p->ionE / U_H;
		double potFrac52 = pow(potentialFrac, 5.0/2.0);
		double potFrac32 = pow(potentialFrac, 3.0/2.0);
		double wQST, eField, eFieldRatio;

		p->rho[0] = electrons;
		for (int i = 0; i < num_t - 1; i++)
		{
			eField = real(p->ee_p[i] + p->ee_m[i]);
			eFieldRatio = abs(eField / E_a);
			if (eFieldRatio != 0) {
				wQST = 4.0 * potFrac52 * nu_a / eFieldRatio * exp(-2.0/3.0 * potFrac32 / eFieldRatio);
			}
			else {
				wQST = 0.0;
			}
			
			change = ht * wQST * neutrals;
			electrons += change;
			neutrals -= change;
			// Don't allow neutrals dip below zero
			if (neutrals < 0.0) neutrals = 0.0;
			if (electrons > num_atoms) electrons = num_atoms;

			p->rho[i + 1] = electrons;


			/* if (eField != 0) {
				p->j_e[i + 1] = p->ionE * wQST * neutrals / eField;
			}
			else {
				p->j_e[i + 1] = 0.0;
			} */
		}

		p->j_e[0] = j_e0;
		fv1 = j_e0;
		for (int i = 0; i < num_t - 1; i++)
		{
			//p->j_e[i + 1] = (1.0 - ht / tauCollision)*p->j_e[i] + ht * pow(charge_e, 2) / mass_e * p->rho[i] * real(p->eFieldPlus[i] + p->eFieldMinus[i]);
			// Andrew original
			//p->j_e[i + 1] = p->j_e[i] * exp(-1.0 * ht / tauCollision) + pow(charge_e, 2) / mass_e * 0.5 * ht * (p->rho[i] * real(p->ee_p[i] + p->ee_m[i]) * exp(-1.0 * ht / tauCollision) + p->rho[i + 1] * real(p->ee_p[i + 1] + p->ee_m[i + 1]));
			
			
			// CALCULATING THE CURRENT for UPPE (based on Ewan's notes Jan-24 2018)
			//  First order method	
			//	current_change = delta_t*(pow(cnst_e,2)/cnst_me)*electrons*real(aptr[t])-delta_t*ve*current;	
			// Second order method (Kolja)
			fv2 = real(p->ee_p[i+1] + p->ee_m[i+1]);
			current_change = ht * (pow(charge_e, 2) / mass_e) * electrons * (fv1 + fv2) * 0.5 - ht * ve * real(p->j_e[i]);
			fv1 = fv2;
			p->j_e[i + 1] = p->j_e[i] + current_change;
		}
		

		//applyWindow(p->j_e);
		fftw_execute(p->np_f);
		normalizeFFT(nl_p, 1);
	}
	else if (p->doPlasmaCalc == 0){
		#pragma omp parallel for
		for (int i = 0; i < num_t; i++)
		{
			p->nl_p[i] = 0.0;
		}

	}
	else {
		cout << "ERROR: Invalid plasma parameter passed to func." << endl;
		exit(EXIT_FAILURE);
	}
	#endif
		
	for (int i = 0; i < freqLowerCutoff; i++) {
		f[i] = 0.0;
		f[i + num_tOver2 + 1] = 0.0;
		f[i + num_t + 2] = 0.0;
		f[i + 3 * num_tOver2 + 3] = 0.0;
	}
	for (int i = freqUpperCutoff + 1; i <= num_tOver2; i++) {
		f[i] = 0.0;
		f[i + num_tOver2 + 1] = 0.0;
		f[i + num_t + 2] = 0.0;
		f[i + 3 * num_tOver2 + 3] = 0.0;
	}
	#pragma omp parallel for
	for (int i = freqLowerCutoff; i <= freqUpperCutoff; i++) {
		//complex<double> deltazA = (1.0i*pow(p->omega[i], 2) / (2.0*(p->k[i])*clightSquared)*p->nl_k[i] + p->omega[i] / (2.0*(p->k[i])*clightSquared*epsilon_0)*p->nl_p[i])*exp(-1.0i*real(p->k[i])*z)*exp(-1.0*abs(imag(p->k[i]))*z);
		complex<double> deltazA = (1.0i*pow(p->omega[i], 2) / (2.0*(p->k[i])*clightSquared)*p->nl_k[i] + p->omega[i] / (2.0*(p->k[i])*clightSquared*epsilon_0)*p->nl_p[i])*exp(-1.0i*(p->k[i])*z);
		f[i] = real(deltazA);
		f[i + num_tOver2 + 1] = imag(deltazA);
		//deltazA = -(1.0i*pow(p->omega[i], 2) / (2.0*(p->k[i])*clightSquared)*p->nl_k[i] + p->omega[i] / (2.0*(p->k[i])*clightSquared*epsilon_0)*p->nl_p[i])*exp(1.0i*real(p->k[i])*z)*exp(-1.0*abs(imag(p->k[i]))*z);
		deltazA = -(1.0i*pow(p->omega[i], 2) / (2.0*(p->k[i])*clightSquared)*p->nl_k[i] + p->omega[i] / (2.0*(p->k[i])*clightSquared*epsilon_0)*p->nl_p[i])*exp(1.0i*(p->k[i])*z);
		f[i + num_t + 2] = real(deltazA);
		f[i + 3 * num_tOver2 + 3] = imag(deltazA); 

	}
	

	return GSL_SUCCESS;
}



void write_multicolumnMonitor(int iterationNo, double theZpos, double *y, odeparam_type *p) {

	char buffer[STRING_BUFFER_SIZE];
	snprintf(buffer, sizeof(char) * STRING_BUFFER_SIZE, "%sPointMon_iter_%i_Zpos_%dnm.dat", SIM_DATA_OUTPUT, iterationNo, (int)round(theZpos * 1.0e9));
	if (VERBOSE >= 7) {
		printf("\t\tWriting point monitor data to file: %s \n", buffer);
	}
	
	int num_tOver2 = num_t/2;
	for (int i = 0; i <= num_tOver2; i++)
	{
		const complex<double> phaseFactor = exp(1.0i * real(p->k[i]) * theZpos) * exp(-1.0 * abs(imag(p->k[i])) * theZpos);
		p->ee_p[i] = (y[i] + 1.0i * y[i + num_tOver2 + 1]) * phaseFactor;
		p->ee_m[num_t - i - 1] = (y[i + num_t + 2] - 1.0i * y[i + 3 * num_tOver2 + 3]) * phaseFactor;

		if (i > 0 && i < num_tOver2) {
			const complex<double> phaseFactor2 = exp(-1.0i * real(p->k[i]) * theZpos) * exp(1.0 * abs(imag(p->k[i])) * theZpos);
			p->ee_p[num_t - i] = (y[i] - 1.0i * y[i + num_tOver2 + 1]) * phaseFactor2;
			p->ee_m[i] = (y[i + num_t + 2] + 1.0i * y[i + 3 * num_tOver2 + 3]) * phaseFactor2;
		}
		
	}

	fftw_execute(p->ep_b);
	fftw_execute(p->em_b);
	normalizeFFT(p->ee_p, 2);
	normalizeFFT(p->ee_m, 2);

	FILE* fp;
	//errno_t err;
	fp = fopen(buffer, "w");
	if (fp != NULL)
	{
		fprintf(fp, "# Values at Monitor zPosition %.17g[micron]  \n", theZpos*1e6);
		fprintf(fp, "# time [sec]\t\tReEtP[V/m]\t\tImEtP[V/m]\t\tReEtM[V/m]\t\tImEtM[V/m]\t\tNumberElectrons\t\tCurrent\n");
		double dt = domain_t / double(num_t);
		for (int i = 0; i < num_t; i++)
		{
			fprintf(fp, "%.7e\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\n", i * dt, real(p->ee_p[i]), imag(p->ee_p[i]), real(p->ee_m[i]), imag(p->ee_m[i]), ne[i], real(j_e[i]));
		}
	}
	else {
		printf("The file %s was not opened\n", buffer);
		exit(-1);
	}
	if (fp != NULL) { fclose(fp); }

	return;
}
