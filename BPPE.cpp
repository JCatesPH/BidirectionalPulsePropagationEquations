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
double sampleLayerThickness, I_0, A_0, tau, lambda_0, omega_0;
double twoColorSH_amplitude, twoColorSH_phase;
int num_t, freqLowerCutoff, freqUpperCutoff, numActiveOmega, numActiveOmega2, l_0;
double domain_t, zStepMaterial1, alpha_tukey;
double *window;
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
	ym0_init = (complex<double>*)malloc(sizeof(complex<double>) * numActiveOmega);
	ym1_init = (complex<double>*)malloc(sizeof(complex<double>) * numActiveOmega);
	ym1_temp = (complex<double>*)malloc(sizeof(complex<double>) * numActiveOmega);
	f0 = (complex<double>*)malloc(sizeof(complex<double>) * numActiveOmega);
	f1 = (complex<double>*)malloc(sizeof(complex<double>) * numActiveOmega);
	integral = (complex<double>*)malloc(sizeof(complex<double>) * numActiveOmega);
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
	//DELME_ArgonDispersion(omegaArray);
	//DELME_AndrewPreformed(omegaArray);
	createWindowFunc(alpha_tukey);
	//if (fp != NULL) { fclose(fp);  }
	set_guess(eFieldPlus, yp_init, ym0_init, ym1_init,ym1_temp,y,eFieldPlusForwardFFT,eFieldMinus,eFieldMinusBackwardFFT,eFieldPlusBackwardFFT, integral);
	
	std::string reldatpath = SIM_DATA_OUTPUT;
	myStructure.writeStructureLayoutToASCIIFile(reldatpath + "StructureLayout.txt");
	myStructure.writeStructureToDATFile(reldatpath + "Structure.dat");
	myMaterialsDB.writeMaterialDBToASCIIFile(reldatpath + "MaterialDatabase.txt");
	myStructure.writeBoundaryLayoutToASCIIFile(reldatpath + "BoundaryLayout.txt");

	/// Early termination
	//printf("############ WARNING ###########:  Early Termination\n"); exit(-1);
    printf("\n\n############ INFO ###########:  ENTERING THE NONLINEAR PART ##############################\n");

	doNonlinearPartofBPPE();

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


void doNonlinearPartofBPPE()
{
	param_type* params = fill_params(chi2_Material1, chi3_Material1, omegaArray, kx, 
			ne, j_e, myStructure.m_layers.begin()->getMaterial().getK(), eFieldPlus, eFieldMinus, nl_k, nl_p, 
			nkForwardFFT, eFieldPlusBackwardFFT, eFieldMinusBackwardFFT, npForwardFFT, myStructure.m_layers.front().getMaterial().getdoPlasmaCalc());
	//gsl_odeiv2_system sys = { func, NULL, (size_t)(4 * numActiveOmega), params };
	//gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk4, zStepMaterial1, epsabs, epsrel);
	//gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, hstart, epsabs, epsrel);
	//gsl_odeiv2_driver_set_nmax(d, ode_nmax);

	const gsl_odeiv2_step_type * stepType = gsl_odeiv2_step_rk4;
	gsl_odeiv2_step * gslStep = gsl_odeiv2_step_alloc(stepType, 2 * num_t + 4);
	gsl_odeiv2_control * gslControl = gsl_odeiv2_control_y_new(epsabs, epsrel);
	gsl_odeiv2_evolve * gslEvolve = gsl_odeiv2_evolve_alloc(2 * num_t + 4);
	gsl_odeiv2_system sys = { func, NULL, (size_t)(4 * numActiveOmega), params };

	double nonlinear_time_initial, nonlinear_time;

	int numzReports, numZsteps = 0;
	double zRight, zStepSize;
	vector<double> endPoints;

	for (int Iteration_number = 1; Iteration_number <= num_iterations; Iteration_number++)
	{
		if (VERBOSE >= 3) { cout << endl << "Iteration number = " << Iteration_number << endl; }
		nonlinear_time_initial = omp_get_wtime();
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

				if (VERBOSE >= 6) { cout << endl << " Doing Layer #" << lit->getlayerIDnum() << endl; }
				for (std::size_t atZpos = 0; atZpos < endPoints.size(); ++atZpos)  {
					zRight = endPoints[atZpos];
					while (zPosition < zRight) 
					{

						GSLerrorFlag = gsl_odeiv2_evolve_apply(gslEvolve, gslControl, gslStep, &sys, &zPosition, zRight, &zStepSize, y);
						if (GSLerrorFlag == GSL_SUCCESS) {
							integrate(zPosition, zStepSize, params, y, integral);
							//printf("I = %d, step = %d, z = %.8g\n", Iteration_number, numZsteps, zPosition);
							numZsteps++;
						}
						else {
							printf("error: driver returned %d\n", GSLerrorFlag);
							break;
						}

						nonlinear_time = omp_get_wtime() - nonlinear_time_initial;
						if ((int)(nonlinear_time / 20) > numzReports) {
							printf("  I = %d, step = %d, z = %.8g, t = %d s\n", Iteration_number, numZsteps, zPosition, (int)nonlinear_time);
							numzReports++;
						}
					}

					printf("  Outputting Point Monitor file at Z location %d[nm]... \n", (int)round(zPosition * 1.0e9));
					//write_multicolumnMonitor(Iteration_number, zPosition, eFieldPlus, eFieldMinus, ne, j_e);
					//params->k = myStructure.m_layers.back().getMaterial().getK();
					write_multicolumnMonitor_Jalen(Iteration_number, zPosition, y, params);
					//params->k = lit->getMaterial().getK();
					
				}

				endPoints.clear();

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

		p = new_initial_data(ym0_init, ym1_init, ym1_temp, yp_init, y, integral);

		nonlinear_time = omp_get_wtime() - nonlinear_time_initial;
		printf("Iteration completed in %.2f seconds with %d steps.\n", nonlinear_time, numZsteps);

		if (Iteration_number > 3) VERBOSE--;
	}
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

	double preformedDensity = 4.0e25;
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

	for (int i = 0; i < freqLowerCutoff; i++) {
		y[i] = 0.0;
		y[i + numActiveOmega] = 0.0;
		y[i + 2 * numActiveOmega] = 0.0;
		y[i + 3 * numActiveOmega] = 0.0;
	}
	for (int i = freqUpperCutoff; i < numActiveOmega; i++) {
		y[i] = 0.0;
		y[i + numActiveOmega] = 0.0;
		y[i + 2 * numActiveOmega] = 0.0;
		y[i + 3 * numActiveOmega] = 0.0;
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

void createWindowFunc(double alpha){
	// Implements a Tukey window
	for (int i = 0; i < num_t/4; i++) {
		window[i] = 0.0;
	}
	for (int i=num_t/4; i < (int)((1+alpha)*num_t/4); i++) {
		window[i] = 0.5 * (1.0 - cos(4*M_PI*(i-num_t/4)/(alpha*num_t)));
	}
	for (int i=(int)((1+alpha)*num_t/4); i <= num_t/2; i++) {
		window[i] = 1.0;
	}
	for (int i=num_t/4; i <= num_t/2; i++) {
		window[num_t-i] = window[i];
	}
	for (int i = 3*num_t/4; i < num_t; i++) {
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

void set_guess(complex<double>* ee_p, complex<double>* yp_init, complex<double>* ym0_init, complex<double>* ym1_init, complex<double>* ym1_temp, double* y, fftw_plan ep_f, complex<double>* ee_m, fftw_plan em_b, fftw_plan ep_b, complex<double>* integral) {

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

		normalizeFFT(ee_p, 1);

		for (int i = 0; i < numActiveOmega; i++) yp_init[i] = ee_p[i];

		writeInputSpectrum(yp_init);
	}

	initalizeArrays(ym1_init, ym0_init, integral);
	fillYfromYpAndYm(y, yp_init, ym0_init); // was initalizeYarray(y, yp_init, ym0_init) BUT they have the same function body
	std::cout << "In SetGuess FINISHED Initalizing stuff" << endl << endl;

	if (VERBOSE >= 7) { cout << "  Going FORWARD through layers" << endl; }
	myStructure.doForwardPassThroughAllBoundaries(y);
	o = 1;

	am_to_zero(y);
	if (VERBOSE >= 7) { cout << "Setting am_to_zero()" << endl; }

	if (VERBOSE >= 7) { cout << "  Going BACKWARD through layers" << endl; }
	myStructure.doBackwardPassThroughAllBoundaries(y);

	update_guess(yp_init, f0, ym1_init, y, integral);

	for (int bb = 1; bb < 3; bb++)
	{
		if (VERBOSE >= 7) { cout << "bb = " << bb << endl; }
		write_out_eFieldAndSpectrumAtZlocation(0, 0, y, 0.0, ee_m, myStructure.m_layers.front().getMaterial().getK(), em_b);

		if (VERBOSE >= 7) { cout << "  Going FORWARD through layers" << endl; }
		myStructure.doForwardPassThroughAllBoundaries(y);
		o = 1;

		write_out_eFieldAndSpectrumAtZlocation(0, 1, y, myStructure.getThickness(), ee_p, myStructure.m_layers.back().getMaterial().getK(), ep_b);
		am_to_zero(y);
		if (VERBOSE >= 7) { cout << "Setting am_to_zero()" << endl; }

		if (VERBOSE >= 7) { cout << "  Going BACKWARD through layers" << endl; }
		myStructure.doBackwardPassThroughAllBoundaries(y);

		int x;
		x = new_initial_data(ym0_init, ym1_init, ym1_temp, yp_init, y, integral);
		
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
		fprintf(fp, "epsabs       \t%.17g\t/*absolute error tolerance for Runge-Kutta */\n", epsabs);
		fprintf(fp, "epsrel       \t%.17g\t/*relative error tolerance for Runge-Kutta */\n", epsrel);
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
			ee[i] = (y[i] + 1.0i * y[i + num_t / 2 + 1]) * exp(-1.0i * real(k[i]) * z) * exp(-1.0 * abs(imag(k[i])) * z);
		}

		for (int i = 1; i < num_t / 2; i++)
		{
			//ee[num_t - i] = (y[i] - 1.0i*y[i + num_t / 2 + 1])*exp(1.0i*conj(k[i])*(zPosition + RHSbufferLayerThickness));		// phase corrections by Andrew (2021-01-25)
			ee[num_t - i] = (y[i] - 1.0i * y[i + num_t / 2 + 1]) * exp(1.0i * real(k[i]) * z) * exp(-1.0 * abs(imag(k[i])) * z);
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
	for (int i = 0; i < freqLowerCutoff; i++) {
		y[i] = 0.0;
		y[i + numActiveOmega] = 0.0;
		y[i + 2 * numActiveOmega] = 0.0;
		y[i + 3 * numActiveOmega] = 0.0;
	}
	for (int i = freqUpperCutoff; i < numActiveOmega; i++) {
		y[i] = 0.0;
		y[i + numActiveOmega] = 0.0;
		y[i + 2 * numActiveOmega] = 0.0;
		y[i + 3 * numActiveOmega] = 0.0;
	}
	return;
}

int new_initial_data(complex<double>*ym0_init, complex<double>*ym1_init, complex<double>*ym1_temp, complex<double>*yp_init, double*y, complex<double>*integral) {

	int nExcRaised = 0;
	double diffNorm = 0.0;
	double vecNorm = 0.0;
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
				ym1_init[i] = ym0_init[i];
			}
			else{
				ym1_init[i] = ((ym0_init[i] * f1[i] - ym1_init[i] * f0[i]) / (f1[i] - f0[i] - ym1_init[i] + ym0_init[i]));
			}		
		}

		// Calculate the norm of the change in backward-prop
		if (normType == -1) {
			// This is the inf-norm
			if (vecNorm < abs(ym1_init[i])) {
				vecNorm = abs(ym1_init[i]);
			}
			if (diffNorm < abs(ym1_init[i] - ym1_temp[i])) {
				diffNorm = abs(ym1_init[i] - ym1_temp[i]);
			}
		}
		else {
			if (nExcRaised == 0) {
				cout << "Exception 22. The norm specified in BPPE.h has not been implemented." << endl;
				nExcRaised++;
			}
		}
	}

	printf("The %d-norm estimate of the rel. err.: %.3e\n", normType, diffNorm/vecNorm);

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

	for (int i = 0; i < freqLowerCutoff; i++) {
		y[i] = 0.0;
		y[i + numActiveOmega] = 0.0;
		y[i + 2 * numActiveOmega] = 0.0;
		y[i + 3 * numActiveOmega] = 0.0;
	}
	for (int i = freqUpperCutoff; i < numActiveOmega; i++) {
		y[i] = 0.0;
		y[i + numActiveOmega] = 0.0;
		y[i + 2 * numActiveOmega] = 0.0;
		y[i + 3 * numActiveOmega] = 0.0;
	}

	for (int i = 0; i <= num_t / 2; i++)
	{
		f0[i] = f1[i];
		ym0_init[i] = ym1_temp[i];
	}
	return 1;
}



int func(double z, const double y[], double f[], void *params) {

	param_type *p = reinterpret_cast<param_type*>(params);

	const int num_tOver2 = num_t / 2;
	const double clightSquared = pow(cLight, 2);
	const double num_td = (double)num_t;

#pragma omp parallel for
	for (int i = 0; i <= num_tOver2; i++)
	{
		const complex<double> phaseFactor = exp(-1.0i * real(p->k[i]) * z) * exp(-1.0 * abs(imag(p->k[i])) * z);
		p->ee_p[i] = (y[i] + 1.0i * y[i + num_tOver2 + 1]) * phaseFactor;
		p->ee_m[num_t - i] = (y[i + num_t + 2] - 1.0i * y[i + 3 * num_tOver2 + 3]) * phaseFactor;

		/*const complex<double> phaseFactorP = exp(-1.0i * p->k[i] * z);
		const complex<double> phaseFactorM = exp(1.0i * p->k[i] * z);
		p->ee_p[i] = (y[i] + 1.0i * y[i + num_tOver2 + 1]) * phaseFactorP;
		p->ee_m[num_t - i] = (y[i + num_t + 2] - 1.0i * y[i + 3 * num_tOver2 + 3]) * phaseFactorM;*/
		// Jalen thinks this is off by 1 due to remainder in the above for loop dividing by 2 
		//p->ee_m[num_t - i - 1] = (y[i + num_t + 2] - 1.0i * y[i + 3 * num_tOver2 + 3]) * phaseFactor;
		// Orignal Andrew was
		

		if (i > 0 && i < num_tOver2) {
			const complex<double> phaseFactor2 = exp(1.0i * real(p->k[i]) * z) * exp(-1.0 * abs(imag(p->k[i])) * z);
			p->ee_p[num_t - i] = (y[i] - 1.0i * y[i + num_tOver2 + 1]) * phaseFactor2;
			p->ee_m[i] = (y[i + num_t + 2] + 1.0i * y[i + 3 * num_tOver2 + 3]) * phaseFactor2;
			//p->ee_p[num_t - i] = (y[i] - 1.0i * y[i + num_tOver2 + 1]) * phaseFactorP;
			//p->ee_m[i] = (y[i + num_t + 2] + 1.0i * y[i + 3 * num_tOver2 + 3]) * phaseFactorM;
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
	normalizeFFT(p->ee_p, 2);
	normalizeFFT(p->ee_m, 2);

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
		complex<double> deltazA = -(1.0i*pow(p->omega[i], 2) / (2.0*(p->k[i])*clightSquared)*p->nl_k[i] + p->omega[i] / (2.0*(p->k[i])*clightSquared*epsilon_0)*p->nl_p[i])*exp(1.0i*real(p->k[i])*z)*exp(-1.0*abs(imag(p->k[i]))*z);
		f[i] = real(deltazA);
		f[i + num_tOver2 + 1] = imag(deltazA);
		deltazA = (1.0i*pow(p->omega[i], 2) / (2.0*(p->k[i])*clightSquared)*p->nl_k[i] + p->omega[i] / (2.0*(p->k[i])*clightSquared*epsilon_0)*p->nl_p[i])*exp(-1.0i*real(p->k[i])*z)*exp(-1.0*abs(imag(p->k[i]))*z);
		f[i + num_t + 2] = real(deltazA);
		f[i + 3 * num_tOver2 + 3] = imag(deltazA); 

	}
	

	return GSL_SUCCESS;
}


void integrate(double z, double zStep, param_type *pars, double*y, complex<double>*integral) {
	const int num_tOver2 = num_t / 2;
	for (int i = 0; i <= num_tOver2; i++)
	{
		if (i < freqLowerCutoff || i > freqUpperCutoff)
		{
			integral[i] = 0.0;
		}
		else {
			//integral[i] += (-1.0i*pow(pars->omega[i], 2) / (2.0*(pars->k[i])*pow(cLight, 2))*nl_k[i] - pars->omega[i] / (2.0*(pars->k[i])*pow(cLight, 2)*epsilon_0)*pars->nl_p[i])* exp(-1.0i*real(pars->k[i]) * z)*exp(-1.0*abs(imag(pars->k[i]))*z)*zStepMaterial1;
			integral[i] += -(1.0i*pow(pars->omega[i], 2) / (2.0*(pars->k[i])*pow(cLight, 2))*pars->nl_k[i] + pars->omega[i] / (2.0*(pars->k[i])*pow(cLight, 2)*epsilon_0)*pars->nl_p[i]) * exp(-1.0i*real(pars->k[i]) * z)*exp(-1.0*abs(imag(pars->k[i]))*z) * zStep;
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
			fprintf(fp, "%.7e\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\n", i * dt, real(eep[i]), imag(eep[i]), real(eem[i]), imag(eem[i]), ne[i], real(j_e[i]));
		}
	}
	else {
		printf("The file %s was not opened\n", buffer);
		exit(-1);
	}
	if (fp != NULL) { fclose(fp); }

	return;
}

void write_multicolumnMonitor_Jalen(int iterationNo, double theZpos, double *y, param_type *p) {

	char buffer[STRING_BUFFER_SIZE];
	snprintf(buffer, sizeof(char) * STRING_BUFFER_SIZE, "%sPointMon_iter_%i_Zpos_%dnm.dat", SIM_DATA_OUTPUT, iterationNo, (int)round(theZpos * 1.0e9));
	printf("\t\tWriting point monitor data to file: %s \n", buffer);

	int num_tOver2 = num_t/2;
	for (int i = 0; i <= num_tOver2; i++)
	{
		const complex<double> phaseFactor = exp(-1.0i * real(p->k[i]) * theZpos) * exp(-1.0 * abs(imag(p->k[i])) * theZpos);
		p->ee_p[i] = (y[i] + 1.0i * y[i + num_tOver2 + 1]) * phaseFactor;
		p->ee_m[num_t - i - 1] = (y[i + num_t + 2] - 1.0i * y[i + 3 * num_tOver2 + 3]) * phaseFactor;

		if (i > 0 && i < num_tOver2) {
			const complex<double> phaseFactor2 = exp(1.0i * real(p->k[i]) * theZpos) * exp(1.0 * abs(imag(p->k[i])) * theZpos);
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

int func_drude(double z, const double y[], double f[], void *params) {
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
		const complex<double> phaseFactor = 1.0;
		p->ee_p[i] = (y[i] + 1.0i * y[i + num_tOver2 + 1]) * phaseFactor;
		p->ee_m[num_t - i] = (y[i + num_t + 2] - 1.0i * y[i + 3 * num_tOver2 + 3]) * phaseFactor;

		/*const complex<double> phaseFactorP = exp(-1.0i * p->k[i] * z);
		const complex<double> phaseFactorM = exp(1.0i * p->k[i] * z);
		p->ee_p[i] = (y[i] + 1.0i * y[i + num_tOver2 + 1]) * phaseFactorP;
		p->ee_m[num_t - i] = (y[i + num_t + 2] - 1.0i * y[i + 3 * num_tOver2 + 3]) * phaseFactorM;*/
		// Jalen thinks this is off by 1 due to remainder in the above for loop dividing by 2 
		//p->ee_m[num_t - i - 1] = (y[i + num_t + 2] - 1.0i * y[i + 3 * num_tOver2 + 3]) * phaseFactor;
		// Orignal Andrew was
		

		if (i > 0 && i < num_tOver2) {
			//const complex<double> phaseFactor2 = exp(1.0i * real(p->k[i]) * z) * exp(-1.0 * abs(imag(p->k[i])) * z);
			const complex<double> phaseFactor2 = 1.0;
			p->ee_p[num_t - i] = (y[i] - 1.0i * y[i + num_tOver2 + 1]) * phaseFactor2;
			p->ee_m[i] = (y[i + num_t + 2] + 1.0i * y[i + 3 * num_tOver2 + 3]) * phaseFactor2;
			//p->ee_p[num_t - i] = (y[i] - 1.0i * y[i + num_tOver2 + 1]) * phaseFactorP;
			//p->ee_m[i] = (y[i + num_t + 2] + 1.0i * y[i + 3 * num_tOver2 + 3]) * phaseFactorM;
		}
		
	}

	
	double omega_pe2 = rho_0 * pow(charge_e, 2) / (mass_e * epsilon_0);
		
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
		//complex<double> deltazA = 1.0i * omega_pe2 / (2.0 * p->k[i] * clightSquared) * (p->ee_p[i] + p->ee_m[i] * exp(2.0i * p->k[i] * z));
		complex<double> deltazA = -omega_pe2 * p->omega[i] / (2.0 * p->k[i] * clightSquared) / (1.0/tauCollision - 1.0i*p->omega[i]) * (p->ee_p[i] + p->ee_m[i] * exp(2.0i * p->k[i] * z));
		f[i] = real(deltazA);
		f[i + num_tOver2 + 1] = imag(deltazA);
		//deltazA = -1.0i * omega_pe2 / (2.0 * p->k[i] * clightSquared) * (p->ee_p[i] * exp(-2.0i * p->k[i] * z) + p->ee_m[i]);
		deltazA = omega_pe2 * p->omega[i] / (2.0 * p->k[i] * clightSquared) / (1.0/tauCollision - 1.0i*p->omega[i]) * (p->ee_p[i] * exp(-2.0i * p->k[i] * z) + p->ee_m[i]);
		f[i + num_t + 2] = real(deltazA);
		f[i + 3 * num_tOver2 + 3] = imag(deltazA); 
		/*
		f[i] = real((-1.0i*pow(p->omega[i], 2) / (2.0*(p->k[i])*clightSquared)*p->nl_k[i] - p->omega[i] / (2.0*(p->k[i])*clightSquared*epsilon_0)*p->nl_p[i])*exp(1.0i*p->k[i]*z));
		f[i + num_tOver2 + 1] = imag((-1.0i * pow(p->omega[i], 2) / (2.0 * (p->k[i]) * clightSquared) * p->nl_k[i] - p->omega[i] / (2.0 * (p->k[i]) * clightSquared * epsilon_0) * p->nl_p[i]) * exp(1.0i * p->k[i] * z));
		f[i + num_t + 2] = real((1.0i * pow(p->omega[i], 2) / (2.0 * (p->k[i]) * clightSquared) * p->nl_k[i] + p->omega[i] / (2.0 * (p->k[i]) * clightSquared * epsilon_0) * p->nl_p[i]) * exp(-1.0i * p->k[i] * z));
		f[i + 3 * num_tOver2 + 3] = imag((1.0i * pow(p->omega[i], 2) / (2.0 * (p->k[i]) * clightSquared) * p->nl_k[i] + p->omega[i] / (2.0 * (p->k[i]) * clightSquared * epsilon_0) * p->nl_p[i]) * exp(-1.0i * p->k[i] * z));
		*/
	}
	

	return GSL_SUCCESS;
}
