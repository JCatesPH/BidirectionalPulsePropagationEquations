#include "BPPE.h"

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

	fftw_execute(p->ep_b);
	fftw_execute(p->em_b);
	normalizeFFT(p->ee_p);
	normalizeFFT(p->ee_m);

	#pragma omp parallel for
	for (int i = 0; i < num_t; i++)
	{
		//p->ee_p[i] = p->ee_p[i] / num_td;
		//p->ee_m[i] = p->ee_m[i] / num_td;
		p->nl_k[i] = p->chi_2 * pow(real(p->ee_p[i] + p->ee_m[i]), 2) + p->chi_3 * pow(real(p->ee_p[i] + p->ee_m[i]), 3);
	}
	
	fftw_execute(p->nk_f);
	normalizeFFT(p->nl_k);

	fftw_execute(p->ep_f);
	fftw_execute(p->em_f);
	normalizeFFT(p->ee_p);
	normalizeFFT(p->ee_m);

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
		//p->nl_k[i] = (-omeg_p2 / (pow(p->omega[i], 2) + 1.0i * p->omega[i] / tauCollision)) * (p->ee_p[i] + p->ee_m[i]);
		//p->nl_k[i] = (-omeg_p2 * tauCollision / (pow(p->omega[i], 2) * tauCollision + 1.0i * p->omega[i])) * (p->ee_p[i] + p->ee_m[i]);
		//p->nl_k[i] = p->nl_k[i] / sqrt((double)num_t);
		//p->nl_k[i] = p->nl_k[i] / sqrt(2.0*M_PI);
		//p->nl_k[i] = ( -omeg_p2 / (p->omega[i] * (p->omega[i] + 1.0i / tauCollision)) ) * (p->ee_p[i] + p->ee_m[i]);
		p->nl_k[i] += ( -omeg_p2 / (p->omega[i] * (p->omega[i] + 1.0i / tauCollision)) ) * (p->ee_p[i] + p->ee_m[i]);
	}

	//fftw_execute(p->nk_f);
	//normalizeFFT(nl_k, 1);

	// ---------------------------------
	#else

	fftw_execute(p->ep_b);
	fftw_execute(p->em_b);
	normalizeFFT(p->ee_p);
	normalizeFFT(p->ee_m);

	#pragma omp parallel for
	for (int i = 0; i < num_t; i++)
	{
		//p->ee_p[i] = p->ee_p[i] / num_td;
		//p->ee_m[i] = p->ee_m[i] / num_td;
		p->nl_k[i] = p->chi_2 * pow(real(p->ee_p[i] + p->ee_m[i]), 2) + p->chi_3 * pow(real(p->ee_p[i] + p->ee_m[i]), 3);
	}
	
	fftw_execute(p->nk_f);
	normalizeFFT(p->nl_k);

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
		normalizeFFT(p->nl_p);

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
		normalizeFFT(p->nl_p);
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



void integrate(double z, double zStep, odeparam_type *pars, double*y, complex<double>*integral) {
	for (int i = 0; i < numActiveOmega; i++)
	{
		if (i < freqLowerCutoff || i > freqUpperCutoff)
		{
			integral[i] = 0.0;
		}
		else {
			integral[i] += -(1.0i*pow(pars->omega[i], 2) / (2.0*(pars->k[i])*pow(cLight, 2))*pars->nl_k[i] + pars->omega[i] / (2.0*(pars->k[i])*pow(cLight, 2)*epsilon_0)*pars->nl_p[i])* exp(-1.0i*real(pars->k[i]) * z)*exp(-1.0*abs(imag(pars->k[i]))*z)*zStep;
		}
	}
	return;
}
