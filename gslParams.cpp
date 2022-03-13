#include "gslParams.h"

ODEParams::ODEParams(int nT, double *omg) { // Constructor definition
	// Takes parameters to set values that are global.
	numT = nT;
	numOmeg = numT / 2 + 1;
	omega = omg;

	// Allocates necessary double arrays
	rho = (double*)malloc(sizeof(double)*numT);
	y = (double*)malloc(sizeof(double)*(4*numOmeg));

	// Allocate complex arrays
	ee_p = (complex<double>*)malloc(sizeof(complex<double>)*numT);
	ee_m = (complex<double>*)malloc(sizeof(complex<double>)*numT);
	nl_k = (complex<double>*)malloc(sizeof(complex<double>)*numT);
	nl_p = (complex<double>*)malloc(sizeof(complex<double>)*numT);
	j_e = (complex<double>*)malloc(sizeof(complex<double>)*numT);

	// Allocate fftw plans
	nk_f = fftw_plan_dft_1d(numT, reinterpret_cast<fftw_complex*>(&nl_k[0]), reinterpret_cast<fftw_complex*>(&nl_k[0]), FFTW_FORWARD, FFTW_WISDOM_TYPE );
	np_f = fftw_plan_dft_1d(numT, reinterpret_cast<fftw_complex*>(&j_e[0]), reinterpret_cast<fftw_complex*>(&nl_p[0]), FFTW_FORWARD, FFTW_WISDOM_TYPE );
	ep_b = fftw_plan_dft_1d(numT, reinterpret_cast<fftw_complex*>(&ee_p[0]), reinterpret_cast<fftw_complex*>(&ee_p[0]), FFTW_BACKWARD, FFTW_WISDOM_TYPE );
	em_b = fftw_plan_dft_1d(numT, reinterpret_cast<fftw_complex*>(&ee_m[0]), reinterpret_cast<fftw_complex*>(&ee_m[0]), FFTW_BACKWARD, FFTW_WISDOM_TYPE );
	ep_f = fftw_plan_dft_1d(numT, reinterpret_cast<fftw_complex*>(&ee_p[0]), reinterpret_cast<fftw_complex*>(&ee_p[0]), FFTW_FORWARD, FFTW_WISDOM_TYPE );
	em_f = fftw_plan_dft_1d(numT, reinterpret_cast<fftw_complex*>(&ee_m[0]), reinterpret_cast<fftw_complex*>(&ee_m[0]), FFTW_FORWARD, FFTW_WISDOM_TYPE );


	// Initialize arrays
	for (int i = 0; i < numT; i++) {
		rho[i] = 0.0;

		ee_p[i] = 0.0;
		ee_m[i] = 0.0;
		nl_k[i] = 0.0;
		nl_p[i] = 0.0;
		j_e[i] = 0.0;
	}

	for (int i = 0; i < numOmeg; i++) {
		y[i] = 0.0;
		y[i + numOmeg] = 0.0;
		y[i + 2*numOmeg] = 0.0;
		y[i + 3*numOmeg] = 0.0;
	}
}

void ODEParams::fillParams(Material myMat) {
	// Load plasma parameters
	doPlasmaCalc = myMat.getdoPlasmaCalc();
	mpi_sigmaK = myMat.getmpi_sigmaK();
	mpi_k = myMat.getmpi_k();
	ionE = myMat.getIonizationEnergy();
	sigmaBremsstrahlung = myMat.getsigmaBremsstrahlung();
	recombTime = myMat.getrecombTime();

	// Load susceptibilities
	chi_2 = myMat.getChi2();
	chi_3 = myMat.getChi3();

	// Load wave numbers
	k = myMat.getK();
	// PLACEHOLDER -> Add kx for 2D calculations
}

ODEParams::ODEParams(const ODEParams &p1) {  // Copy constructor
	// Copy ints and doubles
	numT = p1.numT;
	numOmeg = p1.numOmeg;
	doPlasmaCalc = p1.doPlasmaCalc;
	chi_2 = p1.chi_2;
	chi_3 = p1.chi_3;
	mpi_sigmaK = p1.mpi_sigmaK;
	mpi_k = p1.mpi_k;
	ionE = p1.ionE;
	
	// Copy pure pointers
	omega = p1.omega;

	// Allocates arrays on heap
	rho = (double*)malloc(sizeof(double)*numT);
	y = (double*)malloc(sizeof(double)*(4*numOmeg));
	ee_p = (complex<double>*)malloc(sizeof(complex<double>)*numT);
	ee_m = (complex<double>*)malloc(sizeof(complex<double>)*numT);
	nl_k = (complex<double>*)malloc(sizeof(complex<double>)*numT);
	nl_p = (complex<double>*)malloc(sizeof(complex<double>)*numT);
	j_e = (complex<double>*)malloc(sizeof(complex<double>)*numT);

	// Allocate fftw plans
	nk_f = fftw_plan_dft_1d(numT, reinterpret_cast<fftw_complex*>(&nl_k[0]), reinterpret_cast<fftw_complex*>(&nl_k[0]), FFTW_FORWARD, FFTW_WISDOM_TYPE );
	np_f = fftw_plan_dft_1d(numT, reinterpret_cast<fftw_complex*>(&j_e[0]), reinterpret_cast<fftw_complex*>(&nl_p[0]), FFTW_FORWARD, FFTW_WISDOM_TYPE );
	ep_b = fftw_plan_dft_1d(numT, reinterpret_cast<fftw_complex*>(&ee_p[0]), reinterpret_cast<fftw_complex*>(&ee_p[0]), FFTW_BACKWARD, FFTW_WISDOM_TYPE );
	em_b = fftw_plan_dft_1d(numT, reinterpret_cast<fftw_complex*>(&ee_m[0]), reinterpret_cast<fftw_complex*>(&ee_m[0]), FFTW_BACKWARD, FFTW_WISDOM_TYPE );
	ep_f = fftw_plan_dft_1d(numT, reinterpret_cast<fftw_complex*>(&ee_p[0]), reinterpret_cast<fftw_complex*>(&ee_p[0]), FFTW_FORWARD, FFTW_WISDOM_TYPE );
	em_f = fftw_plan_dft_1d(numT, reinterpret_cast<fftw_complex*>(&ee_m[0]), reinterpret_cast<fftw_complex*>(&ee_m[0]), FFTW_FORWARD, FFTW_WISDOM_TYPE );

	for (int i = 0; i < numT; i++) {
		rho[i] = p1.rho[i];

		ee_p[i] = p1.ee_p[i];
		ee_m[i] = p1.ee_m[i];
		nl_k[i] = p1.nl_k[i];
		nl_p[i] = p1.nl_p[i];
		j_e[i] = p1.j_e[i];
	}

	for (int i = 0; i < numOmeg; i++) {
		y[i] = p1.y[i];
		y[i + numOmeg] = p1.y[i + numOmeg];
		y[i + 2*numOmeg] = p1.y[i + 2*numOmeg];
		y[i + 3*numOmeg] = p1.y[i + 3*numOmeg];
	}

}


void ODEParams::initializeY(complex<double> *sL) { // Function to initialize the y array with LHS source
	for (int i = 0; i < numOmeg; i++)
	{
		y[i] = real(sL[i]);
		y[i + numOmeg] = imag(sL[i]);
		y[i + 2*numOmeg] = 0.0;
		y[i + 3*numOmeg] = 0.0;
	}
}
