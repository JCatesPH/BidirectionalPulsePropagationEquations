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

//#define FFTW_WISDOM_TYPE FFTW_ESTIMATE
#define FFTW_WISDOM_TYPE FFTW_PATIENT

double dtime = omp_get_wtime();
using namespace std;

// Create GLOBAL structure and Material database
Structure myStructure("myFirstStructure");
Simulation mySimulation("andrewApplication1");
MaterialDB myMaterialsDB("myFirstMaterialDB");

// This block of vars were orignally inside main()
double zPosition = 0.0;
double sampleLayerThickness;
vector<double> monitorZlocations;
int GSLerrorFlag, p, oFlag, VERBOSE;
double *omegaArray, *timeValuesArray, *kx, *ne, *y;
complex<double>*eFieldPlus, *eFieldMinus, *yp_init, *ym0_init, *ym1_init, *ym1_temp, *f0, *f1, *integral, *nl_k, *nl_p, *j_e;
fftw_plan nkForwardFFT, eFieldPlusForwardFFT, eFieldPlusBackwardFFT, eFieldMinusForwardFFT, eFieldMinusBackwardFFT, intBackwardFFT, npForwardFFT;
char *paramFileBuffer, SIM_DATA_OUTPUT[30];

int delmeFLAG = 0;

int main(int argc, char *argv[])
{	
	// Setting parameters using input file.
	paramFileBuffer = readParmetersFileToBuffer(argv[1]);
	VERBOSE = getIntParameterValueByName("Verbosity");
	sampleLayerThickness = getDoubleParameterValueByName("sampleLayerThickness");
	getStringParameterValueByName("outputPath", SIM_DATA_OUTPUT);


	//feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#ifdef _WIN64
#ifndef NDEBUG
	_clearfp();
	_controlfp(_controlfp(0, 0) & ~(_EM_INVALID | _EM_ZERODIVIDE | _EM_OVERFLOW),
		_MCW_EM);
#endif
#endif 
	time_t now = time(0);
	char *datetime = ctime(&now);

	if (VERBOSE >= 0) {
		cout << "BPPE code - ACMS Ver.2" << endl;
		cout << "The current date and time: " << datetime << endl;
		cout << "Verbosity = "<< VERBOSE << endl << endl;
		cout << "Working Directory = " << get_current_dir() << SIM_DATA_OUTPUT << endl << endl;
	}

	omp_set_num_threads(num_Threads);
	cout << "Num threads set to  = " << num_Threads << endl << "  Test: ";

	#pragma omp parallel 
	{
		cout << omp_get_thread_num() << " ";
	}
	cout << endl << endl;

	writeSimParameters();


	generateLayers(myMaterialsDB, myStructure);
    setupPointMonitorLocations(myMaterialsDB, myStructure);
	/// VERY Early termination
	//printf("!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!: VERY Early Termination\n"); exit(-1);

	printf("Structure generated and point monitor locations set..\n");
	
	yp_init = (complex<double>*)malloc(sizeof(complex<double>) * numActiveOmega);
	ym0_init = (complex<double>*)malloc(sizeof(complex<double>) * numActiveOmega);
	ym1_init = (complex<double>*)malloc(sizeof(complex<double>) * numActiveOmega);
	ym1_temp = (complex<double>*)malloc(sizeof(complex<double>) * numActiveOmega);
	f0 = (complex<double>*)malloc(sizeof(complex<double>) * numActiveOmega);
	f1 = (complex<double>*)malloc(sizeof(complex<double>) * numActiveOmega);
	integral = (complex<double>*)malloc(sizeof(complex<double>) * numActiveOmega);
	y = (double*)malloc(sizeof(double) * 4 * numActiveOmega);


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
		intBackwardFFT = fftw_plan_dft_2d(num_x, num_t, reinterpret_cast<fftw_complex*>(&integral[0]), reinterpret_cast<fftw_complex*>(&integral[0]), FFTW_BACKWARD, FFTW_WISDOM_TYPE );
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
		intBackwardFFT = fftw_plan_dft_1d(num_t, reinterpret_cast<fftw_complex*>(&integral[0]), reinterpret_cast<fftw_complex*>(&integral[0]), FFTW_BACKWARD, FFTW_WISDOM_TYPE );
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
	DELME_ArgonDispersion(omegaArray);
	//DELME_AndrewPreformed(omegaArray);
	//if (fp != NULL) { fclose(fp);  }
	set_guess(eFieldPlus, yp_init, ym0_init, ym1_init,ym1_temp,f0,f1,y,eFieldPlusForwardFFT,eFieldMinus,eFieldMinusBackwardFFT,eFieldPlusBackwardFFT, integral);
	
	std::string reldatpath = SIM_DATA_OUTPUT;
	myStructure.writeStructureLayoutToASCIIFile(reldatpath + "StructureLayout.txt");
	myStructure.writeStructureToDATFile(reldatpath + "Structure.dat");
	myMaterialsDB.writeMaterialDBToASCIIFile(reldatpath + "MaterialDatabase.txt");
	myStructure.writeBoundaryLayoutToASCIIFile(reldatpath + "BoundaryLayout.txt");

	/// Early termination
	//printf("############ WARNING ###########:  Early Termination\n"); exit(-1);
    printf("\n\n############ INFO ###########:  ENTERING THE NONLINEAR PART ##############################\n");

	doNonlinearPartofBPPE();
	//DELME_doNonlinearPartofBPPE_1Layer(); // Testing single layer simulation

	dtime = omp_get_wtime() - dtime;
	if(VERBOSE >= 0) { cout << "Time in seconds is " << dtime << endl; }
	//cout << "Num threads set to  = " << omp_get_num_threads() << endl << endl;
	cout << endl << "Exiting program.." << endl << endl;

    return 0;
}

void setupPointMonitorLocations(MaterialDB& theMaterialDB, Structure& theStructure)
{
	
	//monitorZlocations.push_back(3e-6); // 3 microns in plasma
	//monitorZlocations.push_back(10e-6); // 10 microns in plasma
	//monitorZlocations.push_back(15e-6); // 10 microns in plasma
	//monitorZlocations.push_back(LHSsourceLayerThickness + 100e-6); // 100 microns in plasma

	//monitorZlocations.push_back(theStructure.getThickness() * 0.25);
	monitorZlocations.push_back(theStructure.getThickness() * 0.5);
	//monitorZlocations.push_back(theStructure.getThickness() * 0.75);

}


void doNonlinearPartofBPPE()
{
	param_type* params = fill_params(chi2_Material1, chi3_Material1, omegaArray, kx, 
			ne, j_e, myStructure.m_layers.begin()->getMaterial().getK(), eFieldPlus, eFieldMinus, nl_k, nl_p, 
			nkForwardFFT, eFieldPlusBackwardFFT, eFieldMinusBackwardFFT, npForwardFFT, myStructure.m_layers.front().getMaterial().getdoPlasmaCalc());
	gsl_odeiv2_system sys = { func, NULL, (size_t)(4 * numActiveOmega), params };
	//gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk4, zStepMaterial1, epsabs, epsrel);
	gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, hstart, epsabs, epsrel);
	gsl_odeiv2_driver_set_nmax(d, ode_nmax);

	/*const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rk4;
	gsl_odeiv2_step * v = gsl_odeiv2_step_alloc(T, 2 * num_t + 4);
	gsl_odeiv2_control * mm = gsl_odeiv2_control_y_new(10.0, 10.0);
	gsl_odeiv2_evolve * f = gsl_odeiv2_evolve_alloc(2 * num_t + 4);
	gsl_odeiv2_system sys = { func, NULL, 2 * num_t + 4, params };*/

	double zRight;
	double zStepSize;

	for (int Iteration_number = 1; Iteration_number <= num_iterations; Iteration_number++)
	{
		if (VERBOSE >= 3) { cout << endl << "Iteration number = " << Iteration_number << endl; }
		delmeFLAG = Iteration_number;
		write_out_eFieldAndSpectrumAtZlocation(Iteration_number, 0, y, zPosition, eFieldMinus, myStructure.m_layers.front().getMaterial().getK(), eFieldMinusBackwardFFT);

		cout << "  Going FORWARD through layers" << endl;
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
				zPosition = lit->getStartZpos();

				if (VERBOSE >= 6) { cout << endl << " Doing Layer #" << lit->getlayerIDnum() << " in " << lit->getNumStepsInLayer() << " z Steps" << endl; }
				for (int aZstep = 0; aZstep < lit->getNumStepsInLayer(); aZstep++)
				{
					// Make sure step size resets don't cause layers to grow
					if (zPosition + lit->getStepSize() > lit->getEndZpos()) {
						zRight = lit->getEndZpos();
					}
					else {
						zRight = zPosition + lit->getStepSize();
					}
					zStepSize = zRight - zPosition;
					
					integrate(zPosition, zStepSize, params, y, integral);
					//GSLerrorFlag = gsl_odeiv2_driver_apply_fixed_step(d, &zPosition, lit->getStepSize(), 1, y);
					GSLerrorFlag = gsl_odeiv2_driver_apply(d, &zPosition, zRight, y);
					
					if (VERBOSE >= 6) {
						printf("step # = %d,  z = %.5g\n", aZstep, zPosition);
					}

					if (GSLerrorFlag == GSL_EMAXITER) {
						printf("######## ERROR #######: GSL driver returned GSL_EMAXITER. The maximum number of steps is set to %d.\n", ode_nmax);
						printf("	Occurred at z = %e\n", zPosition);
						// Reset the evolve and step objects.
						gsl_odeiv2_driver_reset(d);
						// Reset the step size to initial guess.
						gsl_odeiv2_driver_reset_hstart(d, hstart);

						// Decrement aZstep so the reset doesn't cause problems. MAY NOT BE NECESSARY!
						aZstep--;
					} 
					else if (GSLerrorFlag != GSL_SUCCESS)
					{
						printf("error: driver returned %d\n", GSLerrorFlag);
						break;
					}
					// PUT MONITOR OUTPUTS HERE
					if (Iteration_number == num_iterations)
					{
						for (std::size_t atZpos = 0; atZpos < monitorZlocations.size(); ++atZpos) {
							if (zPosition >= monitorZlocations[atZpos] && zPosition < monitorZlocations[atZpos] + myStructure.getZstepSizeAtZpos(monitorZlocations[atZpos]))
							{
								printf("Outputting Point Monitor file at Z location %d[nm]... \n", (int)round(zPosition * 1.0e9));
								write_multicolumnMonitor(Iteration_number, zPosition, eFieldPlus, eFieldMinus, ne, j_e);
							}
						}
					}
				}
				integrate(zPosition, zStepSize, params, y, integral);
				if (VERBOSE >= 3) {
					cout << "Layer #" << lit->getlayerIDnum() << " ending z position: " << zPosition << endl;
				}
				
			}
			//  Finally do the lowside boundary of the last assumed Vacuum layer
			if (lit->getHiSideBoundary() == NULL) {
				boundary(lit->getLowSideBoundary()->m_zPos, lit->getLowSideBoundary()->lowSideLayer->getMaterial().getK(), lit->getLowSideBoundary()->hiSideLayer->getMaterial().getK(), y);
			}
		}

		// Write out 
		//if (Iteration_number == num_iterations)
		    write_out_eFieldAndSpectrumAtZlocation(Iteration_number, 1, y, myStructure.getThickness(), eFieldPlus, myStructure.m_layers.back().getMaterial().getK(), eFieldPlusBackwardFFT);
		if (VERBOSE >= 4) { cout << "Setting am_to_zero()" << endl; }
		am_to_zero(y);

		cout << "  Going BACKWARD through BOUNDARIES" << endl;
		myStructure.doBackwardPassThroughAllBoundaries(y);

		p = new_initial_data(ym0_init, ym1_init, ym1_temp, yp_init, f0, f1, y, integral);

		if (Iteration_number > 3) VERBOSE--;
	}
}

void DELME_doNonlinearPartofBPPE_1Layer()
{
	param_type* params = fill_params(chi2_Material1, chi3_Material1, omegaArray, kx, 
			ne, j_e, myStructure.m_layers.begin()->getMaterial().getK(), eFieldPlus, eFieldMinus, nl_k, nl_p, 
			nkForwardFFT, eFieldPlusBackwardFFT, eFieldMinusBackwardFFT, npForwardFFT, myStructure.m_layers.begin()->getMaterial().getdoPlasmaCalc());
	gsl_odeiv2_system sys = { func, NULL, (size_t)(4 * numActiveOmega), params };
	//gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk4, zStepMaterial1, epsabs, epsrel);
	gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, hstart, epsabs, epsrel);
	gsl_odeiv2_driver_set_nmax(d, ode_nmax);

	/*const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rk4;
	gsl_odeiv2_step * v = gsl_odeiv2_step_alloc(T, 2 * num_t + 4);
	gsl_odeiv2_control * mm = gsl_odeiv2_control_y_new(10.0, 10.0);
	gsl_odeiv2_evolve * f = gsl_odeiv2_evolve_alloc(2 * num_t + 4);
	gsl_odeiv2_system sys = { func, NULL, 2 * num_t + 4, params };*/

	double zRight;
	double zStepSize;

	for (int Iteration_number = 1; Iteration_number <= num_iterations; Iteration_number++)
	{
		if (VERBOSE >= 3) { cout << endl << "Iteration number = " << Iteration_number << endl; }
		delmeFLAG = Iteration_number;
		write_out_eFieldAndSpectrumAtZlocation(Iteration_number, 0, y, zPosition, eFieldMinus, myStructure.m_layers.begin()->getMaterial().getK(), eFieldMinusBackwardFFT);

		cout << "  Going FORWARD through layers" << endl;
		// Skip the LHS layer and the RHS layers
		params->k = myStructure.m_layers.begin()->getMaterial().getK();
		params->chi_2 = myStructure.m_layers.begin()->getMaterial().getChi2();
		params->chi_3 = myStructure.m_layers.begin()->getMaterial().getChi3();
		params->doPlasmaCalc = myStructure.m_layers.begin()->getMaterial().getdoPlasmaCalc();
		params->mpi_k = myStructure.m_layers.begin()->getMaterial().getmpi_k();
		params->mpi_sigmaK = myStructure.m_layers.begin()->getMaterial().getmpi_sigmaK();
		zPosition = myStructure.m_layers.begin()->getStartZpos();

		if (VERBOSE >= 6) { cout << endl << " Doing Layer #" << myStructure.m_layers.begin()->getlayerIDnum() << " in " << myStructure.m_layers.begin()->getNumStepsInLayer() << " z Steps" << endl; }
		for (int aZstep = 0; aZstep < myStructure.m_layers.begin()->getNumStepsInLayer(); aZstep++)
	{
			// Make sure step size resets don't cause layers to grow
			if (zPosition + myStructure.m_layers.begin()->getStepSize() > myStructure.m_layers.begin()->getEndZpos()) {
				zRight = myStructure.m_layers.begin()->getEndZpos();
	}
	else {
				zRight = zPosition + myStructure.m_layers.begin()->getStepSize();
	}
			zStepSize = zRight - zPosition;

			integrate(zPosition, zStepSize, params, y, integral);
			//GSLerrorFlag = gsl_odeiv2_driver_apply_fixed_step(d, &zPosition, lit->getStepSize(), 1, y);
			GSLerrorFlag = gsl_odeiv2_driver_apply(d, &zPosition, zRight, y);
			
			if (VERBOSE >= 6) {
				printf("step # = %d,  z = %.5g\n", aZstep, zPosition);
}

			if (GSLerrorFlag == GSL_EMAXITER) {
				printf("######## ERROR #######: GSL driver returned GSL_EMAXITER. The maximum number of steps is set to %d.\n", ode_nmax);
				printf("	Occurred at z = %e\n", zPosition);
				// Reset the evolve and step objects.
				gsl_odeiv2_driver_reset(d);
				// Reset the step size to initial guess.
				gsl_odeiv2_driver_reset_hstart(d, hstart);

				// Decrement aZstep so the reset doesn't cause problems. MAY NOT BE NECESSARY!
				aZstep--;
			} 
			else if (GSLerrorFlag != GSL_SUCCESS)
	{
				printf("error: driver returned %d\n", GSLerrorFlag);
				break;
	}
			// PUT MONITOR OUTPUTS HERE
			if (Iteration_number == num_iterations)
	{
				for (std::size_t atZpos = 0; atZpos < monitorZlocations.size(); ++atZpos) {
					if (zPosition >= monitorZlocations[atZpos] && zPosition < monitorZlocations[atZpos] + myStructure.getZstepSizeAtZpos(monitorZlocations[atZpos]))
					{
						printf("Outputting Point Monitor file at Z location %d[nm]... \n", (int)round(zPosition * 1.0e9));
						write_multicolumnMonitor(Iteration_number, zPosition, eFieldPlus, eFieldMinus, ne, j_e);
	}
				}
			}
		}
		integrate(zPosition, zStepSize, params, y, integral);

		// Write out 
		//if (Iteration_number == num_iterations)
		write_out_eFieldAndSpectrumAtZlocation(Iteration_number, 1, y, myStructure.getThickness(), eFieldPlus, myStructure.m_layers.back().getMaterial().getK(), eFieldPlusBackwardFFT);
		if (VERBOSE >= 4) { cout << "Setting am_to_zero()" << endl; }
		am_to_zero(y);

		cout << "  Going BACKWARD through BOUNDARIES" << endl;
		//myStructure.doBackwardPassThroughAllBoundaries(y);

		p = new_initial_data(ym0_init, ym1_init, ym1_temp, yp_init, f0, f1, y, integral);

		if (Iteration_number > 3) VERBOSE--;
	}
}

void fill_omg_k(double*omg, double*kx) {

/*	if (numDimensionsMinusOne == 1) {
	
		for (int j = 0; j <= num_x / 2; j++)
		{
			if (j == 0) {
				for (int i = 0; i <= num_t / 2; i++)
				{
					kx[i] = (M_PI / domain_x)*j;
					omg[i] = (M_PI / domain_t)*i;
				}
			}
			else if (j == num_x / 2) {
				for (int i = 0; i <= num_t / 2; i++)
				{
					kx[i + numActiveOmega2] = (M_PI / domain_x)*j;
					omg[i + numActiveOmega2] = (M_PI / domain_t)*i;
				}
			}
			else {
				for (int i = 0; i < num_t; i++)
				{
					kx[i + (j - 1)*num_t + (num_t / 2 + 1)] = (M_PI / domain_x)*j;
					if (i <= num_t / 2) {
						omg[i + (j - 1)*num_t + (num_t / 2 + 1)] = (M_PI / domain_t)* (double)i;
					}
					else {
						omg[i + (j - 1)*num_t + (num_t / 2 + 1)] = (M_PI / domain_t)*((double)i - num_t);
					}

				}
			}
		}
		for (int j = 0; j < numActiveOmega; j++)
		{
			if (j == 0) {
				k_0[j] = 0.0;
				k_1[j] = 0.0;
				k_2[j] = 0.0;
				k_3[j] = 0.0;
			}
			else 
			{
				k_0[j] = copysign(1.0, omg[j])*sqrt(pow(omg[j], 2) * pow(index_0(omg[j], kx[j]), 2) / (pow(cLight, 2)) - pow(kx[j], 2));
				k_1[j] = copysign(1.0, omg[j])*sqrt(pow(omg[j], 2) * pow(index_1(omg[j], kx[j], fp), 2) / (pow(cLight, 2)) - pow(kx[j], 2));
				k_2[j] = copysign(1.0, omg[j])*sqrt(pow(omg[j], 2) * pow(index_2(omg[j], kx[j]), 2) / (pow(cLight, 2)) - pow(kx[j], 2));
				k_3[j] = copysign(1.0, omg[j])*sqrt(pow(omg[j], 2) * pow(index_3(omg[j], kx[j]), 2) / (pow(cLight, 2)) - pow(kx[j], 2));
			}		
		} 
	} */
	//else {
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
		/* ANDREW ORIGNAL
		for (int i = 0; i <= num_t / 2; i++)
		{
			if (i == 0) {
				k_0[i] = 0.0;
				k_1[i] = 0.0;
				k_2[i] = 0.0;
				k_3[i] = 0.0;
			}
			else {
				k_0[i] = omg[i] * index_0(omg[i], 0.0) / cLight;
				k_1[i] = omg[i] * index_1(omg[i], 0.0, err, fp) / cLight;
				k_2[i] = omg[i] * index_2(omg[i], 0.0) / cLight;
				k_3[i] = omg[i] * index_3(omg[i], 0.0) / cLight;
			}
		}*/
	//}
	return;
}

#define ANDREW_PREFORMED_METHOD1
void DELME_AndrewPreformed(double* omg) {
	Material *plasmaMat;
	plasmaMat = myMaterialsDB.getMaterialByName("PlasmaMat");

	double preformedDensity = 4.3e25;
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
		plasmaMat->m_k[i] = omg[i] / cLight * sqrt(complex<double>(1.0 - pow(omega_plasma / omg[i], 2)));
		fprintf(fp, "%.7g\t%.17g\t%.17g \n", omg[i], real(plasmaMat->m_k[i]), imag(plasmaMat->m_k[i]));
	}
	if (fp != NULL) { fclose(fp); }
}

void DELME_ArgonDispersion(double* omg) {
	Material *ArgonMat;
	complex<double> n0;
	double lambda;
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
	for (int i = 0; i < numActiveOmega; i++)
	{
		lambda = 2 * M_PI*cLight / omg[i];
		//n0 = 1+2.50141e-3/(91.012-pow(lambda,-2))+5.00283e-4/(87.892-pow(lambda,-2))+5.22343e-2/(214.02-pow(lambda,-2));
		n0 = 1+6.432135E-5+2.8606021E-2/(144-pow(lambda,-2)); // From Peck and Fisher (1964), valid in range: 0.47-2.06 um

		n0 += -1.0i * 1e-7 * omg[i] / (1+omg[i]);
		if (i == 0) {
			ArgonMat->m_k[i] = 0.0;
		}
		else {
			ArgonMat->m_k[i] = omg[i] * n0 / cLight;
		}

		fprintf(fp, "%.7g\t%.17g\t%.17g \n", lambda, real(n0), imag(n0));
	}
	if (fp != NULL) { fclose(fp); }
}

//void copy_omg_k_ToMaterial(MaterialDB aMaterialDB, double* kx, complex<double>* k_0, complex<double>* k_1, complex<double>* k_2, complex<double>* k_3) {
//
//}

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

void writeInputSpectrum(std::complex<double>* yp_init)
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
			fprintf(fp_spectrum, "%g \t %+.17g \t %+.17g \t %+.17g\n", (2 * M_PI / domain_t) * i, real(yp_init[i]), imag(yp_init[i]), abs(yp_init[i]));
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

void initalizeArrays(std::complex<double>* ym1_init, std::complex<double>* ym0_init, std::complex<double>* integral)
{
	for (int i = 0; i < numActiveOmega; i++)
	{
		if (i == 0) {
			ym1_init[i] = 0.0;
			ym0_init[i] = 0.0;
			integral[i] = 0.0;
		}
		else {
			ym1_init[i] = INITIAL_GUESS_SEED_VALUE; //exact_linear(k_0, k_1, k_2, k_0, k_0, i)*yp_init[i];
			ym0_init[i] = 0.0;  //-1.0*yp_init[i] * exp(-1.0i*k_0[i] * 2.0*LHSsourceLayerThickness);
			integral[i] = 0.0;
		}
	}
}

void initalizeYarray(double* y, std::complex<double>* yp_init, std::complex<double>* ym0_init)
{
	for (int i = 0; i < 4 * numActiveOmega; i++)
	{
		if (i < numActiveOmega) {
			y[i] = real(yp_init[i]);
		}
		else if (i < 2 * numActiveOmega) {
			y[i] = imag(yp_init[i - numActiveOmega]);
		}
		else if (i < 3 * numActiveOmega) {
			y[i] = real(ym0_init[i - 2 * numActiveOmega]);
		}
		else {
			y[i] = imag(ym0_init[i - 3 * numActiveOmega]);
		}
	}
}

void fillYfromYpAndYm(double* y, std::complex<double>* yp_init, std::complex<double>* ym0_init)
{
	for (int i = 0; i < 4 * numActiveOmega; i++)
	{
		if (i < numActiveOmega) {
			y[i] = real(yp_init[i]);
		}
		else if (i < 2 * numActiveOmega) {
			y[i] = imag(yp_init[i - numActiveOmega]);
		}
		else if (i < 3 * numActiveOmega) {
			y[i] = real(ym0_init[i - 2 * numActiveOmega]);
		}
		else {
			y[i] = imag(ym0_init[i - 3 * numActiveOmega]);
		}
	}
}

void write2DtoFile(std::complex<double>* ee_p)
{
	char buffer[STRING_BUFFER_SIZE];
	snprintf(buffer, sizeof(char) * STRING_BUFFER_SIZE, "%sEp_2d_%i.dat", SIM_DATA_OUTPUT, 0);
	FILE* fp;
	//errno_t err;
	fp = fopen(buffer, "w");
	if (fp != NULL)
	{
		for (int i = 0; i < num_t * num_x; i++)
		{
			fprintf(fp, "%.10lf \t", real(ee_p[i]));
		}
	}
	else {
		printf("The file 'Ep_2d_0.dat' was not opened\n");
		exit(-1);
	}
}


void set_guess(complex<double>* ee_p, complex<double>* yp_init, complex<double>* ym0_init, complex<double>* ym1_init, complex<double>* ym1_temp, complex<double>* f0, complex<double>* f1, double* y, fftw_plan ep_f, complex<double>* ee_m, fftw_plan em_b, fftw_plan ep_b, complex<double>* integral) {

	double ht = (1.0 * domain_t) / double(num_t);
	int o;
	if (numDimensionsMinusOne == 1) {
		// 2D source code deleted temporarly see orignal set_guess() above
		std::cout << "WARNING WARNING ###### ABORTING #####: In SetGuess 2D code. Source code was deleted temporarly see orignal set_guess()" << endl << endl;
		exit(-1);
	}
	else {
		for (int i = 0; i < num_t; i++) {
			double s = (-domain_t/2.0) + ht * ((double)i + 1);
			ee_p[i] = A_0 * (sqrt(1.0 - twoColorSH_amplitude) * exp(-2.0 * log(2.0) * pow(s, 2) / (pow(tau, 2))) * cos(omega_0 * s)
				           + sqrt(twoColorSH_amplitude) *       exp(-8.0 * log(2.0) * pow(s, 2) / (pow(tau, 2))) * cos(2.0 * omega_0 * s + twoColorSH_phase));
		}
		writeInputEfield(ee_p);
		fftw_execute(ep_f);

		for (int i = 0; i < numActiveOmega; i++) yp_init[i] = ee_p[i];

		writeInputSpectrum(yp_init);
	}

	initalizeArrays(ym1_init, ym0_init, integral);
	fillYfromYpAndYm(y, yp_init, ym0_init); // was initalizeYarray(y, yp_init, ym0_init) BUT they have the same function body
	std::cout << "In SetGuess FINISHED Initalizing stuff" << endl << endl;

#ifdef	USE_CPP_BOUNDARY
	if (VERBOSE >= 7) { cout << "  Going FORWARD through layers" << endl; }
	myStructure.doForwardPassThroughAllBoundaries(y);
	o = 1;
#else
	cout << "  Going FORWARD through layers USING ANDREW?" << endl;
	//boundary(LHSsourceLayerThickness, k_0, k_1, y);
	//for (int aLayer = 0; aLayer <= numLayersInSample - 1; aLayer++)
	// Andrew orignal deleted stuff was here
#endif

	am_to_zero(y);
	if (VERBOSE >= 7) { cout << "Setting am_to_zero()" << endl; }

#ifdef	USE_CPP_BOUNDARY
	if (VERBOSE >= 7) { cout << "  Going BACKWARD through layers" << endl; }
	myStructure.doBackwardPassThroughAllBoundaries(y);
#else
	cout << "  Going BACKWARD through layers USING ANDREW?" << endl;
	//if (o == 1)
	//{
	//	boundary(LHSsourceLayerThickness + lengthSample, k_3, k_2, y);
	//	for (int aLayer = 0; aLayer < numLayersInSample - 1; aLayer++)
	//		// Andrew orignal deleted stuff was here
	//}
	//else {
	//	boundary(LHSsourceLayerThickness + lengthSample, k_3, k_1, y);
	//	for (int aLayer = 0; aLayer < numLayersInSample - 1; aLayer++)
	// 	   	// Andrew orignal deleted stuff was here
	//}
	//boundary(LHSsourceLayerThickness, k_1, k_0, y);
#endif

	update_guess(yp_init, f0, ym1_init, y, integral);

	for (int bb = 1; bb < 3; bb++)
	{
		if (VERBOSE >= 7) { cout << "bb = " << bb << endl; }
		write_out_eFieldAndSpectrumAtZlocation(0, 0, y, 0.0, ee_m, myStructure.m_layers.front().getMaterial().getK(), em_b);

#ifdef	USE_CPP_BOUNDARY
		if (VERBOSE >= 7) { cout << "  Going FORWARD through layers" << endl; }
		myStructure.doForwardPassThroughAllBoundaries(y);
		o = 1;
#else
		//cout << "  Going FOREWARD through layers USING ANDREW??? o=" << o << endl;
		//boundary(LHSsourceLayerThickness, k_0, k_1, y);
		//for (int aLayer = 0; aLayer <= numLayersInSample - 1; aLayer++)
			// Andrew orignal deleted stuff was here
#endif

		write_out_eFieldAndSpectrumAtZlocation(0, 1, y, myStructure.getThickness(), ee_p, myStructure.m_layers.back().getMaterial().getK(), ep_b);
		am_to_zero(y);
		if (VERBOSE >= 7) { cout << "Setting am_to_zero()" << endl; }

#ifdef	USE_CPP_BOUNDARY
		if (VERBOSE >= 7) { cout << "  Going BACKWARD through layers" << endl; }
		myStructure.doBackwardPassThroughAllBoundaries(y);
#else
		//cout << "  Going BACKWARD through layers USING ANDREW??? o=" << o << endl;
		//if (o == 1)
		//{

		//	boundary(LHSsourceLayerThickness + lengthSample, k_3, k_2, y);
		//	for (int aLayer = 0; aLayer < numLayersInSample - 1; aLayer++)
		//	// Andrew orignal deleted stuff was here
		else {
			//boundary(LHSsourceLayerThickness + lengthSample, k_3, k_1, y);
			//for (int aLayer = 0; aLayer < numLayersInSample - 1; aLayer++)
			// 	// Andrew orignal deleted stuff was here
		//}
		//boundary(LHSsourceLayerThickness, k_1, k_0, y);
#endif

		int x;
		if (myStructure.m_layers.size() > 1) {
		x = new_initial_data(ym0_init, ym1_init, ym1_temp, yp_init, f0, f1, y, integral);
	}
		else {
			x = DELME_new_initial_data_1Layer(ym0_init, ym1_init, ym1_temp, yp_init, f0, f1, y, integral);
		}
		
	}

	for (int i = 0; i < numActiveOmega; i++)
	{
		if (i == 0) { ym0_init[i] = 0.0; }
		else { ym0_init[i] = INITIAL_GUESS_SEED_VALUE;	}
	}

	fillYfromYpAndYm(y, yp_init, ym0_init);

	return;
}


param_type* fill_params(double chi_2, double chi_3, double* omg, double* kx, double* ne, complex<double>* j_e, complex<double>* k, complex<double>* ee_p, complex<double>* ee_m, complex<double>* nl_k, complex<double>* nl_p, fftw_plan nk_f, fftw_plan ep_b, fftw_plan em_b, fftw_plan np_f, int plasmaBool) {

	param_type* r = (param_type*)malloc(sizeof(param_type));
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
	for (int i = 1; i <= num_t / 2; i++)
	{
		if (i >= freqLowerCutoff && i <= freqUpperCutoff) {
			double k_0_r = real(k_0[i]);
			double k_0_i = imag(k_0[i]);
			double k_1_r = real(k_1[i]);
			double k_1_i = imag(k_1[i]);

			sp = (exp(-1.0i*((k_0_r + 1.0i*k_0_i) - (k_1_r + 1.0i*k_1_i))*z)*((k_0_r + 1.0i*k_0_i) + (k_1_r + 1.0i*k_1_i)) / (2.0*(k_1_r + 1.0i*k_1_i))*(y[i] + 1.0i*y[i + num_t / 2 + 1]) + exp(1.0i*((k_0_r + 1.0i*k_0_i) + (k_1_r + 1.0i*k_1_i))*z)*((k_1_r + 1.0i*k_1_i) - (k_0_r + 1.0i*k_0_i)) / (2.0*(k_1_r + 1.0i*k_1_i))*(y[i + num_t + 2] + 1.0i*y[i + 3 * num_t / 2 + 3]));

			sm = (exp(-1.0i*((k_0_r + 1.0i*k_0_i) + (k_1_r + 1.0i*k_1_i))*z)*((k_1_r + 1.0i*k_1_i) - (k_0_r + 1.0i*k_0_i)) / (2.0*(k_1_r + 1.0i*k_1_i))*(y[i] + 1.0i*y[i + num_t / 2 + 1]) + exp(1.0i*((k_0_r + 1.0i*k_0_i) - (k_1_r + 1.0i*k_1_i))*z)*((k_0_r + 1.0i*k_0_i) + (k_1_r + 1.0i*k_1_i)) / (2.0*(k_1_r + 1.0i*k_1_i))*(y[i + num_t + 2] + 1.0i*y[i + 3 * num_t / 2 + 3]));

			y[i] = real(sp);
			y[i + num_t / 2 + 1] = imag(sp);
			y[i + num_t + 2] = real(sm);
			y[i + 3 * num_t / 2 + 3] = imag(sm);
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
		fprintf(fp, "freqUpperCutoff \t%.17g\t/* upper frequency cut-off */\n", (double)freqUpperCutoff);
		fprintf(fp, "freqLowerCutoff \t%.17g\t/*lower frequency cut-off */\n", (double)freqLowerCutoff);
		fprintf(fp, "num_iterations  \t%.17g\t/*number of BPPE iterations*/ */\n", (double)num_iterations);
		fprintf(fp, "shift           \t%.17g\t/*FILL */\n", shift);
		fprintf(fp, "numDimensionsMinusOne \t%.17g\t/*(1+1) dimension (0) or (2+1) dimension (1) */\n", (double)numDimensionsMinusOne);
		fprintf(fp, "plasmaOnOff     \t%d\t/*plasma on (1) or off (0) */\n", plasmaOnOff);
		fprintf(fp, "l_0               \t%.17g\t/*FILL */\n", (double)l_0);
		fprintf(fp, "numActiveOmega    \t%.17g\t/*FILL */\n", (double)numActiveOmega);
		fprintf(fp, "numActiveOmega2   \t%.17g\t/*FILL */\n", (double)numActiveOmega2);
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
		fprintf(fp, "Keldysh                \t%.17g\t/*keldysh parameter */\n", Keldysh);
		fprintf(fp, "omegaPlasma            \t%.17g\t/*fully-ionized plasma frequency */\n", omegaPlasma);
		fprintf(fp, "Sellmeir_chi_1_1       \t%.17g\t/*FILL */\n", Sellmeir_chi_1_1);
		fprintf(fp, "Sellmeir_chi_1_2       \t%.17g\t/*FILL */\n", Sellmeir_chi_1_2);
		fprintf(fp, "Sellmeir_chi_1_3       \t%.17g\t/*FILL */\n", Sellmeir_chi_1_3);
		fprintf(fp, "Sellmeir_omega_1       \t%.17g\t/*FILL */\n", Sellmeir_omega_1);
		fprintf(fp, "Sellmeir_omega_2       \t%.17g\t/*FILL */\n", Sellmeir_omega_2);
		fprintf(fp, "Sellmeir_omega_3       \t%.17g\t/*FILL */\n", Sellmeir_omega_3);

		// Print the GSL ODE parameters
		fprintf(fp, "epsabs       \t%.17g\t/*absolute error tolerance for Runge-Kutta */\n", epsabs);
		fprintf(fp, "epsrel       \t%.17g\t/*relative error tolerance for Runge-Kutta */\n", epsrel);
		fprintf(fp, "ode_nmax       \t%.2g\t/*maximum ode steps before reset */\n", (double)ode_nmax);
		fprintf(fp, "hstart      \t%.17g\t/*initial guess of step size*/\n", hstart);
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
			ee[i] = (y[i] + 1.0i * y[i + num_t / 2 + 1]) * exp(-1.0i * k[i] * z);
		}

		for (int i = 1; i < num_t / 2; i++)
		{
			//ee[num_t - i] = (y[i] - 1.0i*y[i + num_t / 2 + 1])*exp(1.0i*conj(k[i])*(zPosition + RHSbufferLayerThickness));		// phase corrections by Andrew (2021-01-25)
			ee[num_t - i] = (y[i] - 1.0i * y[i + num_t / 2 + 1]) * exp(1.0i * conj(k[i]) * z);
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
			if (i < num_t / 2) {
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

	for (int i = 0; i < num_t; i++)
	{
		ee[i] = ee[i] / (double(num_t));
	}

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

// PARIS Sets second half to zero (minuns going real and imag parts)
void am_to_zero(double*y) {

	for (int i = 0; i <= num_t + 1; i++)
	{
		y[i + num_t + 2] = 0.0;
	}

	return;
}

void update_guess(complex<double>*yp_init, complex<double>*f0, complex<double>*ym1_init, double*y, complex<double>*integral) {

	for (int i = 0; i <= num_t / 2; i++)
	{
		f0[i] = (y[i + num_t + 2] + 1.0i*y[i + 3 * num_t / 2 + 3]) + integral[i];
		integral[i] = 0.0;
	}

	for (int i = 0; i < 2 * num_t + 4; i++)
	{
		if (i <= num_t / 2) {
			y[i] = real(yp_init[i]);
		}
		else if (i <= num_t + 1) {
			y[i] = imag(yp_init[i - (num_t / 2 + 1)]);
		}
		else if (i <= 3 * num_t / 2 + 2) {
			y[i] = real(ym1_init[i - (num_t + 2)]);
		}
		else {
			y[i] = imag(ym1_init[i - (3 * num_t / 2 + 3)]);
		}

	}
	return;
}

int new_initial_data(complex<double>*ym0_init, complex<double>*ym1_init, complex<double>*ym1_temp, complex<double>*yp_init, complex<double>*f0, complex<double>*f1, double*y, complex<double>*integral) {

	int nExcRaised = 0;
	for (int i = 0; i <= num_t / 2; i++)
	{
		f1[i] = (y[i + num_t + 2] + 1.0i*y[i + 3 * num_t / 2 + 3]) + integral[i];
		integral[i] = 0.0;
		ym1_temp[i] = ym1_init[i];
	}

	for (int i = 0; i <= num_t / 2; i++)
	{
		if (i == 0) {
			ym1_init[i] = 0.0;
		}
		else {
			// Added exception handling for case where field is not updated. 
			if (abs(f1[i] - f0[i] - ym1_init[i] + ym0_init[i]) < DBL_MIN) {
				if (nExcRaised == 0) {
					cout << "Exception 21. Underflow of denominator in secant method. Further iterations may not be necessary!" << endl;
					nExcRaised++;
				}
				ym1_init[i] = ym0_init[i];
			}
			else{
				ym1_init[i] = ((ym0_init[i] * f1[i] - ym1_init[i] * f0[i]) / (f1[i] - f0[i] - ym1_init[i] + ym0_init[i]));
			}		
		}
	}

	for (int i = 0; i < 2 * num_t + 4; i++)
	{
		if (i <= num_t / 2) {
			y[i] = real(yp_init[i]);
		}
		else if (i <= num_t + 1) {
			y[i] = imag(yp_init[i - (num_t / 2 + 1)]);
		}
		else if (i <= 3 * num_t / 2 + 2) {
			y[i] = real(ym1_init[i - (num_t + 2)]);
		}
		else {
			y[i] = imag(ym1_init[i - (3 * num_t / 2 + 3)]);
		}
	}

	for (int i = 0; i <= num_t / 2; i++)
	{
		f0[i] = f1[i];
		ym0_init[i] = ym1_temp[i];
	}
	return 1;
}

int DELME_new_initial_data_1Layer(complex<double>*ym0_init, complex<double>*ym1_init, complex<double>*ym1_temp, complex<double>*yp_init, complex<double>*f0, complex<double>*f1, double*y, complex<double>*integral) {

	for (int i = 0; i <= num_t / 2; i++)
	{
		f1[i] = (y[i + num_t + 2] + 1.0i*y[i + 3 * num_t / 2 + 3]) + integral[i];
		integral[i] = 0.0;
		ym1_temp[i] = ym1_init[i];
		ym1_init[i] = 0.0;
	}

	for (int i = 0; i < 2 * num_t + 4; i++)
	{
		if (i <= num_t / 2) {
			y[i] = real(yp_init[i]);
		}
		else if (i <= num_t + 1) {
			y[i] = imag(yp_init[i - (num_t / 2 + 1)]);
		}
		else if (i <= 3 * num_t / 2 + 2) {
			y[i] = real(ym1_init[i - (num_t + 2)]);
		}
		else {
			y[i] = imag(ym1_init[i - (3 * num_t / 2 + 3)]);
		}
	}

	for (int i = 0; i <= num_t / 2; i++)
	{
		f0[i] = f1[i];
		ym0_init[i] = ym1_temp[i];
	}
	return 1;
}


int func_correct(double z, const double y[], double f[], void *params) {

	param_type *p = reinterpret_cast<param_type*>(params);

	const int num_tOver2 = num_t / 2;
	const double clightSquared = pow(cLight, 2);
	const double num_td = (double)num_t;

#pragma omp parallel for
	for (int i = 0; i <= num_tOver2; i++)
	{
		//const complex<double> phaseFactor = exp(-1.0i * real(p->k[i]) * z) * exp(-1.0 * abs(imag(p->k[i])) * z);
		const complex<double> phaseFactorP = exp(-1.0i * p->k[i] * z);
		const complex<double> phaseFactorM = exp(1.0i * p->k[i] * z);
		p->ee_p[i] = (y[i] + 1.0i * y[i + num_tOver2 + 1]) * phaseFactorP;
		// Jalen thinks this is off by 1 due to remainder in the above for loop dividing by 2 
		//p->ee_m[num_t - i - 1] = (y[i + num_t + 2] - 1.0i * y[i + 3 * num_tOver2 + 3]) * phaseFactor;
		// Orignal Andrew was
		p->ee_m[num_t - i] = (y[i + num_t + 2] - 1.0i * y[i + 3 * num_tOver2 + 3]) * phaseFactorM;

		if (i > 0 && i < num_tOver2) {
			p->ee_p[num_t - i] = (y[i] - 1.0i * y[i + num_tOver2 + 1]) * phaseFactorP;
			p->ee_m[i] = (y[i + num_t + 2] + 1.0i * y[i + 3 * num_tOver2 + 3]) * phaseFactorM;
		}
		
	}

/* #pragma omp  parallel for
	for (int i = 1; i < num_tOver2; i++)
	{
		const complex<double> phaseFactor2 = exp(1.0i * real(p->k[i]) * z) * exp(-1.0 * abs(imag(p->k[i])) * z);
		p->ee_p[num_t - i] = (y[i] - 1.0i * y[i + num_tOver2 + 1]) * phaseFactor2;
		p->ee_m[i] = (y[i + num_t + 2] + 1.0i * y[i + 3 * num_tOver2 + 3]) * phaseFactor2;
	} */


	fftw_execute(p->ep_b);
	fftw_execute(p->em_b);

#pragma omp parallel for
	for (int i = 0; i < num_t; i++)
	{
		p->ee_p[i] = p->ee_p[i] / num_td;
		p->ee_m[i] = p->ee_m[i] / num_td;
		p->nl_k[i] = p->chi_2 * pow(real(p->ee_p[i] + p->ee_m[i]), 2) + p->chi_3 * pow(real(p->ee_p[i] + p->ee_m[i]), 3);
	}


	fftw_execute(p->nk_f);
	
	if (p->doPlasmaCalc == 2) {

		// POSSIBLE ERROR WHY FACTOR 2.0 in following ht calculation???
		double ht = domain_t / num_td;
		double neutrals = num_atoms - rho_0;                          // Neutral particles
		double electrons = rho_0;                      // background Electrons
		double change = 0.0;    
		double current = 0.0;                                       // Current to be exported to UPPE
		double current_change = 0.0e0;                              // Current change for differential equation
		double ve = 0.0e0, fv1 = 0.0e0, fv2 = 0.0e0;

		p->rho[0] = electrons;
		for (int i = 0; i < num_t - 1; i++)
		{
			double fieldIntensity = (( pow(real(p->ee_p[i] + p->ee_m[i]),2) ) / Znaught );  // MIRO real+real
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
			fv2 = real(p->ee_p[i+1] + p->ee_m[i+1]);
			current_change = ht * (pow(charge_e, 2) / mass_e) * electrons * (fv1 + fv2) * 0.5 - ht * ve * current;
			fv1 = fv2;
			current += current_change;
			p->j_e[i + 1] = current;
		}

		//if (z == distanceSourceToSample)
		

		// I THINK.. this takes the Current j_e and forewardFFTs it into the array called nl_p which is then used in the Integrate()
		fftw_execute(p->np_f);

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
		/* f[i] = real((-1.0i*pow(p->omega[i], 2) / (2.0*(p->k[i])*clightSquared)*p->nl_k[i] - p->omega[i] / (2.0*(p->k[i])*clightSquared*epsilon_0)*p->nl_p[i])*exp(1.0i*real(p->k[i])*z)*exp(-1.0*abs(imag(p->k[i]))*z));
		f[i + num_tOver2 + 1] = imag((-1.0i * pow(p->omega[i], 2) / (2.0 * (p->k[i]) * clightSquared) * p->nl_k[i] - p->omega[i] / (2.0 * (p->k[i]) * clightSquared * epsilon_0) * p->nl_p[i]) * exp(1.0i * real(p->k[i]) * z) * exp(-1.0 * abs(imag(p->k[i])) * z));
		f[i + num_t + 2] = real((1.0i * pow(p->omega[i], 2) / (2.0 * (p->k[i]) * clightSquared) * p->nl_k[i] + p->omega[i] / (2.0 * (p->k[i]) * clightSquared * epsilon_0) * p->nl_p[i]) * exp(-1.0i * real(p->k[i]) * z) * exp(-1.0 * abs(imag(p->k[i])) * z));
		f[i + 3 * num_tOver2 + 3] = imag((1.0i * pow(p->omega[i], 2) / (2.0 * (p->k[i]) * clightSquared) * p->nl_k[i] + p->omega[i] / (2.0 * (p->k[i]) * clightSquared * epsilon_0) * p->nl_p[i]) * exp(-1.0i * real(p->k[i]) * z) * exp(-1.0 * abs(imag(p->k[i])) * z)); */
		f[i] = real((-1.0i*pow(p->omega[i], 2) / (2.0*(p->k[i])*clightSquared)*p->nl_k[i] - p->omega[i] / (2.0*(p->k[i])*clightSquared*epsilon_0)*p->nl_p[i])*exp(1.0i*p->k[i]*z));
		f[i + num_tOver2 + 1] = imag((-1.0i * pow(p->omega[i], 2) / (2.0 * (p->k[i]) * clightSquared) * p->nl_k[i] - p->omega[i] / (2.0 * (p->k[i]) * clightSquared * epsilon_0) * p->nl_p[i]) * exp(1.0i * p->k[i] * z));
		f[i + num_t + 2] = real((1.0i * pow(p->omega[i], 2) / (2.0 * (p->k[i]) * clightSquared) * p->nl_k[i] + p->omega[i] / (2.0 * (p->k[i]) * clightSquared * epsilon_0) * p->nl_p[i]) * exp(-1.0i * p->k[i] * z));
		f[i + 3 * num_tOver2 + 3] = imag((1.0i * pow(p->omega[i], 2) / (2.0 * (p->k[i]) * clightSquared) * p->nl_k[i] + p->omega[i] / (2.0 * (p->k[i]) * clightSquared * epsilon_0) * p->nl_p[i]) * exp(-1.0i * p->k[i] * z));
	}
	

	return GSL_SUCCESS;
}


void integrate(double z, double zStep, param_type *pars, double*y, complex<double>*integral) {
	const int num_tOver2 = num_t / 2;

/* #pragma omp parallel for	
	for (int i = 0; i <= num_tOver2; i++)
	{
		pars->ee_p[i] = (y[i] + 1.0i*y[i + num_tOver2 + 1])*exp(-1.0i*real(pars->k[i])*z)*exp(-1.0*abs(imag(pars->k[i]))*z);
	}

#pragma omp parallel for	
	for (int i = 1; i < num_tOver2; i++)
	{
		pars->ee_p[num_t - i] = (y[i] - 1.0i*y[i + num_tOver2 + 1])*exp(1.0i*real(pars->k[i])*z)*exp(-1.0*abs(imag(pars->k[i]))*z);
	}

#pragma omp parallel for	
	for (int i = 0; i <= num_tOver2; i++)
	{
		pars->ee_m[i] = (y[i + num_t + 2] + 1.0i*y[i + 3 * num_tOver2 + 3])*exp(1.0i*real(pars->k[i])*z)*exp(-1.0*abs(imag(pars->k[i]))*z);
	}

#pragma omp parallel for	
	for (int i = 1; i < num_tOver2; i++)
	{
		pars->ee_m[num_t - i] = (y[i + num_t + 2] - 1.0i*y[i + 3 * num_tOver2 + 3])*exp(-1.0i*real(pars->k[i])*z)*exp(-1.0*abs(imag(pars->k[i]))*z);
	}

	fftw_execute(pars->ep_b);
	fftw_execute(pars->em_b);

#pragma omp parallel for	
	for (int i = 0; i < num_t; i++)
	{
		pars->ee_p[i] = pars->ee_p[i] / (double(num_t));
		pars->ee_m[i] = pars->ee_m[i] / (double(num_t));
		pars->nl_k[i] = pars->chi_2 * pow(real(pars->ee_p[i] + pars->ee_m[i]), 2) + pars->chi_3 * pow(real(pars->ee_p[i] + pars->ee_m[i]), 3);
		}

	fftw_execute(pars->nk_f);
	
	if (pars->doPlasmaCalc == 2) {
		double ht = (2.0 * domain_t) / double(num_t);
		double neutrals = num_atoms;                          // Neutral particles
		double electrons = 0.0; 
		double change = 0.0;
		double current = 0.0;                                       // Current to be exported to UPPE
		double current_change = 0.0e0;                              // Current change for differential equation
		double ve = 0.0e0, fv1 = 0.0e0, fv2 = 0.0e0;

		pars->rho[0] = electrons;
		for (int i = 0; i < num_t - 1; i++)
		{
			double fieldIntensity = ((pow(real(pars->ee_p[i] + pars->ee_m[i]), 2)) / Znaught);
			change = neutrals * pars->mpi_sigmaK * ht * pow(fieldIntensity, pars->mpi_k);
			electrons += change;
			neutrals -= change;
			// Don't allow neutrals dip below zero
			if (neutrals < 0.0) neutrals = 0.0;
			if (electrons > num_atoms) electrons = num_atoms;

			pars->rho[i + 1] = electrons;
		}

		j_e[0] = j_e0;
		fv1 = j_e0;
		for (int i = 0; i < num_t - 1; i++)
		{
			// CALCULATING THE CURRENT for UPPE (based on Ewan's notes Jan-24 2018)
			//  First order method	
			//	current_change = delta_t*(pow(cnst_e,2)/cnst_me)*electrons*real(aptr[t])-delta_t*ve*current;	
			// Second order method (Kolja)
			fv2 = real(pars->ee_p[i + 1] + pars->ee_m[i + 1]);
			current_change = ht * (pow(charge_e, 2) / mass_e) * electrons * (fv1 + fv2) * 0.5 - ht * ve * current;
			fv1 = fv2;
			current += current_change;
			j_e[i + 1] = current;
		}

		fftw_execute(pars->np_f);
	}
	else if (pars->doPlasmaCalc == 0) {
		for (int i = 0; i < num_t; i++)
		{
			nl_p[i] = 0.0;
		}
	}
	else {
		cout << "ERROR: Invalid option passed for plasma calculation." << endl;
		exit(EXIT_FAILURE);
	}
 */
	for (int i = 0; i <= num_tOver2; i++)
	{
		if (i < freqLowerCutoff || i > freqUpperCutoff)
		{
			integral[i] = 0.0;
		}
		else {
			//integral[i] += (-1.0i*pow(pars->omega[i], 2) / (2.0*(pars->k[i])*pow(cLight, 2))*nl_k[i] - pars->omega[i] / (2.0*(pars->k[i])*pow(cLight, 2)*epsilon_0)*pars->nl_p[i])* exp(-1.0i*real(pars->k[i]) * z)*exp(-1.0*abs(imag(pars->k[i]))*z)*zStepMaterial1;
			integral[i] += (-1.0i*pow(pars->omega[i], 2) / (2.0*(pars->k[i])*pow(cLight, 2))*nl_k[i] - pars->omega[i] / (2.0*(pars->k[i])*pow(cLight, 2)*epsilon_0)*pars->nl_p[i]) * exp(-1.0i * pars->k[i] * z) * zStep;
		}
	}
	return;
}


void write_multicolumnMonitor(int iterationNo, double theZpos, complex<double>* eep, complex<double>* eem, double* ne, complex<double>* j_e) {

	char buffer[STRING_BUFFER_SIZE];
	snprintf(buffer, sizeof(char) * STRING_BUFFER_SIZE, "%sPointMon_iter_%i_Zpos_%dnm.dat", SIM_DATA_OUTPUT, iterationNo, (int)round(theZpos * 1.0e9));
	printf("\t\tWriting point monitor data to file: %s \n", buffer);

	FILE* fp;
	//errno_t err;
	fp = fopen(buffer, "w");
	if (fp != NULL)
	{
		fprintf(fp, "# Algorithm iteration number %d  at Moitor zPosition %.17g[micron]  \n", iterationNo,  theZpos*1e6);
		fprintf(fp, "# time [sec]\t\tReEtP[V/m]\t\tImEtP[V/m]\t\tReEtM[V/m]\t\tImEtM[V/m]\t\tNumberElectrons\t\tCurrent\n");
		double dt = domain_t / double(num_t);
		for (int i = 0; i < num_t; i++)
		{
			fprintf(fp, "%.7e\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\n", i * dt, real(eep[i]), imag(eep[i]), real(eem[i]), imag(eep[i]), ne[i], real(j_e[i]));
		}
	}
	else {
		printf("The file %s was not opened\n", buffer);
		exit(-1);
	}
	if (fp != NULL) { fclose(fp); }

	return;
}

int func(double z, const double y[], double f[], void *params) {
	/*
	================================
			WARNING!!
	================================
	This version of func is experimental and only meant to reproduce "Method 2" of Andrew's code.
	*/

	param_type *p = reinterpret_cast<param_type*>(params);

	const int num_tOver2 = num_t / 2;
	const double clightSquared = pow(cLight, 2);
	const double num_td = (double)num_t;

#pragma omp parallel for
	for (int i = 0; i <= num_tOver2; i++)
	{
		//const complex<double> phaseFactor = exp(-1.0i * real(p->k[i]) * z) * exp(-1.0 * abs(imag(p->k[i])) * z);
		const complex<double> phaseFactorP = exp(1.0i * p->k[i] * z);
		const complex<double> phaseFactorM = exp(-1.0i * p->k[i] * z);
		p->ee_p[i] = (y[i] + 1.0i * y[i + num_tOver2 + 1]) * phaseFactorM;
		// Jalen thinks this is off by 1 due to remainder in the above for loop dividing by 2 
		//p->ee_m[num_t - i - 1] = (y[i + num_t + 2] - 1.0i * y[i + 3 * num_tOver2 + 3]) * phaseFactor;
		// Orignal Andrew was
		p->ee_m[num_t - i] = (y[i + num_t + 2] - 1.0i * y[i + 3 * num_tOver2 + 3]) * phaseFactorP;

		if (i > 0 && i < num_tOver2) {
			p->ee_p[num_t - i] = (y[i] - 1.0i * y[i + num_tOver2 + 1]) * phaseFactorP;
			p->ee_m[i] = (y[i + num_t + 2] + 1.0i * y[i + 3 * num_tOver2 + 3]) * phaseFactorM;
		}
		
	}

/* #pragma omp  parallel for
	for (int i = 1; i < num_tOver2; i++)
	{
		const complex<double> phaseFactor2 = exp(1.0i * real(p->k[i]) * z) * exp(-1.0 * abs(imag(p->k[i])) * z);
		p->ee_p[num_t - i] = (y[i] - 1.0i * y[i + num_tOver2 + 1]) * phaseFactor2;
		p->ee_m[i] = (y[i + num_t + 2] + 1.0i * y[i + 3 * num_tOver2 + 3]) * phaseFactor2;
	} */


	fftw_execute(p->ep_b);
	fftw_execute(p->em_b);

#pragma omp parallel for
	for (int i = 0; i < num_t; i++)
	{
		p->ee_p[i] = p->ee_p[i] / num_td;
		p->ee_m[i] = p->ee_m[i] / num_td;
		p->nl_k[i] = p->chi_2 * pow(real(p->ee_p[i] + p->ee_m[i]), 2) + p->chi_3 * pow(real(p->ee_p[i] + p->ee_m[i]), 3);
	}


	fftw_execute(p->nk_f);
	
	// POSSIBLE ERROR WHY FACTOR 2.0 in following ht calculation???
	double ht = domain_t / num_td;
	double change, w_QST;
	const double U_i = 15.76; // Ionization potential of Argon [eV]
	const double U_H = 13.6; // Ionization potential of Hydrgogen [eV]
	const double E_a = 5.14e11; // ? [V/m]
	const double mu_a = 4.13e16; // ? [Hz]
	double ve = 0.0e0, fv1 = 0.0e0, fv2 = 0.0e0;
	const double charge2 = pow(charge_e, 2);

	p->rho[0] = rho_0;
	for (int i = 0; i < num_t - 1; i++)
	{
		w_QST = 4 * mu_a * pow(U_i/U_H, 5/2) / abs(p->ee_p[i] / E_a) * exp(-2 * pow(U_i/U_H, 3/2) / (3 * abs(p->ee_p[i]/E_a)));
		change = w_QST * (num_atoms - p->rho[i]);
		p->rho[i + 1] = ht * change + p->rho[i];
	}

	for (int i = 0; i < num_t; i++)
	{
		p->j_e[i] = p->rho[i] * real(p->ee_p[i]);
	}

	fftw_execute(p->np_f);
	for (int i = 0; i < num_t; i++)
	{
		p->nl_p[i] = mu_0 * charge2 / mass_e * p->nl_p[i];
	}

		
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
		/* f[i] = real((-1.0i*pow(p->omega[i], 2) / (2.0*(p->k[i])*clightSquared)*p->nl_k[i] - p->omega[i] / (2.0*(p->k[i])*clightSquared*epsilon_0)*p->nl_p[i])*exp(1.0i*real(p->k[i])*z)*exp(-1.0*abs(imag(p->k[i]))*z));
		f[i + num_tOver2 + 1] = imag((-1.0i * pow(p->omega[i], 2) / (2.0 * (p->k[i]) * clightSquared) * p->nl_k[i] - p->omega[i] / (2.0 * (p->k[i]) * clightSquared * epsilon_0) * p->nl_p[i]) * exp(1.0i * real(p->k[i]) * z) * exp(-1.0 * abs(imag(p->k[i])) * z));
		f[i + num_t + 2] = real((1.0i * pow(p->omega[i], 2) / (2.0 * (p->k[i]) * clightSquared) * p->nl_k[i] + p->omega[i] / (2.0 * (p->k[i]) * clightSquared * epsilon_0) * p->nl_p[i]) * exp(-1.0i * real(p->k[i]) * z) * exp(-1.0 * abs(imag(p->k[i])) * z));
		f[i + 3 * num_tOver2 + 3] = imag((1.0i * pow(p->omega[i], 2) / (2.0 * (p->k[i]) * clightSquared) * p->nl_k[i] + p->omega[i] / (2.0 * (p->k[i]) * clightSquared * epsilon_0) * p->nl_p[i]) * exp(-1.0i * real(p->k[i]) * z) * exp(-1.0 * abs(imag(p->k[i])) * z)); */
		f[i] = real((-1.0i*pow(p->omega[i], 2) / (2.0*(p->k[i])*clightSquared)*p->nl_k[i] - p->omega[i] / (2.0*(p->k[i])*clightSquared*epsilon_0)*p->nl_p[i])*exp(1.0i*p->k[i]*z));
		f[i + num_tOver2 + 1] = imag((-1.0i * pow(p->omega[i], 2) / (2.0 * (p->k[i]) * clightSquared) * p->nl_k[i] - p->omega[i] / (2.0 * (p->k[i]) * clightSquared * epsilon_0) * p->nl_p[i]) * exp(1.0i * p->k[i] * z));
		f[i + num_t + 2] = real((1.0i * pow(p->omega[i], 2) / (2.0 * (p->k[i]) * clightSquared) * p->nl_k[i] + p->omega[i] / (2.0 * (p->k[i]) * clightSquared * epsilon_0) * p->nl_p[i]) * exp(-1.0i * p->k[i] * z));
		f[i + 3 * num_tOver2 + 3] = imag((1.0i * pow(p->omega[i], 2) / (2.0 * (p->k[i]) * clightSquared) * p->nl_k[i] + p->omega[i] / (2.0 * (p->k[i]) * clightSquared * epsilon_0) * p->nl_p[i]) * exp(-1.0i * p->k[i] * z));
	}
	

	return GSL_SUCCESS;
}