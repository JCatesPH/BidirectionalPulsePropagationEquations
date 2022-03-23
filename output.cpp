#include "BPPE.h"

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

		}
	}
	else {
		printf("Failed to open file '%s'\n", inputEfieldFilePathName2);

		exit(-1);
	}
	if (fp_spectrum != NULL) { fclose(fp_spectrum); }
}


void write_out_eFieldAndSpectrumAtZlocation(int num, int j, double*y, double z, complex<double>*ee, complex<double>*k, fftw_plan e_b) {

	char efieldFilePathName[STRING_BUFFER_SIZE];
	char spectrumFilePathName[STRING_BUFFER_SIZE];
	// COLMTODO E_ij -> E_iter_i_Reflect(0) ,Transmitted(1)
	snprintf(efieldFilePathName, sizeof(char) * STRING_BUFFER_SIZE, "%sEfield_iteration_%i_%s.dat", SIM_DATA_OUTPUT, num, (j == 0 ? "Reflected" : "Transmitted"));
	snprintf(spectrumFilePathName, sizeof(char) * STRING_BUFFER_SIZE, "%sSpectrum_iteration_%i_%s.dat", SIM_DATA_OUTPUT, num, (j == 0 ? "Reflected" : "Transmitted"));
	//snprintf(buffer2, sizeof(char) * STRING_BUFFER_SIZE, "Snew_%i%i.dat", Iteration_number, j);

	/* if 	(z>myStructure.getThickness()) {
		printf("Request to write field at point outside structure %g  Structure range 0<->%g\n", z, myStructure.getThickness());
	} */


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
	normalizeFFT(ee, fftnorm);
	//applyWindow(ee);

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



void write_multicolumnMonitor(int iterationNo, double theZpos, double *y, ODEParams *odeObj) {

	char buffer[STRING_BUFFER_SIZE];
	snprintf(buffer, sizeof(char) * STRING_BUFFER_SIZE, "%sPointMon_iter_%i_Zpos_%dnm.dat", SIM_DATA_OUTPUT, iterationNo, (int)round(theZpos * 1.0e9));
	if (VERBOSE >= 7) {
		printf("\t\tWriting point monitor data to file: %s \n", buffer);
	}
	
	int num_tOver2 = num_t/2;
	for (int i = 0; i <= num_tOver2; i++)
	{
		const complex<double> phaseFactor = exp(1.0i * real(odeObj->k[i]) * theZpos) * exp(-1.0 * abs(imag(odeObj->k[i])) * theZpos);
		odeObj->ee_p[i] = (y[i] + 1.0i * y[i + num_tOver2 + 1]) * phaseFactor;
		odeObj->ee_m[num_t - i - 1] = (y[i + num_t + 2] - 1.0i * y[i + 3 * num_tOver2 + 3]) * phaseFactor;

		if (i > 0 && i < num_tOver2) {
			const complex<double> phaseFactor2 = exp(-1.0i * real(odeObj->k[i]) * theZpos) * exp(1.0 * abs(imag(odeObj->k[i])) * theZpos);
			odeObj->ee_p[num_t - i] = (y[i] - 1.0i * y[i + num_tOver2 + 1]) * phaseFactor2;
			odeObj->ee_m[i] = (y[i + num_t + 2] + 1.0i * y[i + 3 * num_tOver2 + 3]) * phaseFactor2;
		}
		
	}

	fftw_execute(odeObj->ep_b);
	fftw_execute(odeObj->em_b);
	normalizeFFT(odeObj->ee_p, fftnorm);
	normalizeFFT(odeObj->ee_m, fftnorm);
	//applyWindow(odeObj->ee_p);
	//applyWindow(odeObj->ee_m);


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
			fprintf(fp, "%.7e\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\n", i * dt, real(odeObj->ee_p[i]), imag(odeObj->ee_p[i]), real(odeObj->ee_m[i]), imag(odeObj->ee_m[i]), odeObj->rho[i], real(odeObj->j_e[i]));
		}
	}
	else {
		printf("The file %s was not opened\n", buffer);
		exit(-1);
	}
	if (fp != NULL) { fclose(fp); }

	return;
}
