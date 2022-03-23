#include "BPPE.h"

int dAdz(double z, const double y[], double f[], void *odep) {

	//odeparam_type *p = reinterpret_cast<odeparam_type*>(odep);
	ODEParams *odeObj = reinterpret_cast<ODEParams*>(odep);

	const int num_tOver2 = num_t / 2;
	const double num_td = (double)num_t;

	#pragma omp parallel for
	for (int i = 0; i < numActiveOmega; i++)
	{
		//const complex<double> phaseFactor = exp(1.0i * real(odeObj->k[i]) * z) * exp(-1.0 * abs(imag(odeObj->k[i])) * z);
		//const complex<double> phaseFactor2 = exp(-1.0i * real(odeObj->k[i]) * z) * exp(-1.0 * abs(imag(odeObj->k[i])) * z);
		const complex<double> phaseFactor = exp(1.0i * odeObj->k[i] * z);
		const complex<double> phaseFactor2 = exp(-1.0i * odeObj->k[i] * z);

		odeObj->ee_p[i] = (y[i] + 1.0i * y[i + numActiveOmega]) * phaseFactor;
		odeObj->ee_m[i] = (y[i + 2 * numActiveOmega] + 1.0i * y[i + 3 * numActiveOmega]) * phaseFactor2;

		if (i > 0 && i < numActiveOmega - 1) {	
			odeObj->ee_p[num_t - i] = (y[i] - 1.0i * y[i + numActiveOmega]) * phaseFactor2;
			odeObj->ee_m[num_t - i] = (y[i + 2 * numActiveOmega] - 1.0i * y[i + 3 * numActiveOmega]) * phaseFactor;	
		}
	}

	// DELETE ME : Drude model
	#ifdef DO_DRUDE_MODEL

	fftw_execute(odeObj->ep_b);
	fftw_execute(odeObj->em_b);
	normalizeFFT(odeObj->ee_p);
	normalizeFFT(odeObj->ee_m);
	applyWindow(odeObj->ee_p);
	applyWindow(odeObj->ee_m);

	#pragma omp parallel for
	for (int i = 0; i < num_t; i++)
	{
		//odeObj->ee_p[i] = odeObj->ee_p[i] / num_td;
		//odeObj->ee_m[i] = odeObj->ee_m[i] / num_td;
		odeObj->p_nl[i] = odeObj->chi_2 * pow(real(odeObj->ee_p[i] + odeObj->ee_m[i]), 2) + odeObj->chi_3 * pow(real(odeObj->ee_p[i] + odeObj->ee_m[i]), 3);
	}
	
	fftw_execute(odeObj->p_ffft);
	normalizeFFT(odeObj->p_nl);

	fftw_execute(odeObj->ep_f);
	fftw_execute(odeObj->em_f);
	normalizeFFT(odeObj->ee_p);
	normalizeFFT(odeObj->ee_m);

	const double sig0 = rho_0 * pow(charge_e, 2) * tauCollision / mass_e;
	const double omeg_p2 = rho_0 * pow(charge_e, 2) / (mass_e * epsilon_0); 
	#pragma omp parallel for
	for (int i = 0; i < num_t; i++)
	{	
		// Calculate the current density
		//odeObj->jhat[i] = sig0 / (1.0 - 1.0i * odeObj->omega[i] * tauCollision) * (odeObj->ee_p[i]+odeObj->ee_m[i]);
		//const double w = sig0 / (1.0 + pow(odeObj->omega[i]*tauCollision, 2));
		//odeObj->jhat[i] = w * (1.0 + 1.0i * odeObj->omega[i] * tauCollision) * (odeObj->ee_p[i]+odeObj->ee_m[i]);
		odeObj->jhat[i] = 0.0;

		// Calculate the polarization
		//odeObj->p_nl[i] = (1.0 - omeg_p2 / (pow(odeObj->omega[i], 2) + 1.0i * odeObj->omega[i] / tauCollision)) * (odeObj->ee_p[i]+odeObj->ee_m[i]);
		//odeObj->p_nl[i] = (-omeg_p2 / (pow(odeObj->omega[i], 2) + 1.0i * odeObj->omega[i] / tauCollision)) * (odeObj->ee_p[i] + odeObj->ee_m[i]);
		//odeObj->p_nl[i] = (-omeg_p2 * tauCollision / (pow(odeObj->omega[i], 2) * tauCollision + 1.0i * odeObj->omega[i])) * (odeObj->ee_p[i] + odeObj->ee_m[i]);
		//odeObj->p_nl[i] = odeObj->p_nl[i] / sqrt((double)num_t);
		//odeObj->p_nl[i] = odeObj->p_nl[i] / sqrt(2.0*M_PI);
		//odeObj->p_nl[i] = ( -omeg_p2 / (odeObj->omega[i] * (odeObj->omega[i] + 1.0i / tauCollision)) ) * (odeObj->ee_p[i] + odeObj->ee_m[i]);
		odeObj->p_nl[i] += ( -omeg_p2 / (odeObj->omega[i] * (odeObj->omega[i] + 1.0i / tauCollision)) ) * (odeObj->ee_p[i] + odeObj->ee_m[i]);
	}

	//fftw_execute(odeObj->p_ffft);
	//normalizeFFT(p_nl, 1);

	// ---------------------------------
	#else

	fftw_execute(odeObj->ep_b);
	fftw_execute(odeObj->em_b);
	normalizeFFT(odeObj->ee_p, fftnorm);
	normalizeFFT(odeObj->ee_m, fftnorm);
	//applyWindow(odeObj->ee_p);
	//applyWindow(odeObj->ee_m);

	#pragma omp parallel for
	for (int i = 0; i < num_t; i++)
	{
		//odeObj->ee_p[i] = odeObj->ee_p[i] / num_td;
		//odeObj->ee_m[i] = odeObj->ee_m[i] / num_td;
		odeObj->p_nl[i] = odeObj->chi_2 * pow(real(odeObj->ee_p[i] + odeObj->ee_m[i]), 2) + odeObj->chi_3 * pow(real(odeObj->ee_p[i] + odeObj->ee_m[i]), 3);
	}
	
	fftw_execute(odeObj->p_ffft);
	normalizeFFT(odeObj->p_nl, fftnorm);

	if (odeObj->doPlasmaCalc == 3) {
		double ht = domain_t / num_td;
		double neutrals = num_atoms - rho_0;      // Neutral particles
		double electrons = rho_0;                 // background Electrons
		double change = 0.0;                      // dN/dt
		double current = 0.0;                     // Current to be exported to UPPE
		double current_change = 0.0e0;            // Current change for differential equation
		double ve = 3.33e14;                      // Electron-ion collision frequency
		double vr = 1 / odeObj->recombTime;       // Recombination frequency
		double gamma = odeObj->sigmaBremsstrahlung / odeObj->ionE;
		double fv1 = 0.0e0, fv2 = 0.0e0;

		odeObj->rho[0] = electrons;
		for (int i = 0; i < num_t - 1; i++)
		{
			double fieldIntensity = ( pow(real(odeObj->ee_p[i] + odeObj->ee_m[i]),2) ) / Znaught ;  // MIRO real+real
			// First, MPI ionization is calculated
			change = neutrals * odeObj->mpi_sigmaK * ht * pow(fieldIntensity, odeObj->mpi_k);
			// Next, the avalanche ionization
			change += ht * gamma * fieldIntensity * electrons;
			// Last, recombination term is computed
			change -= ht * vr * electrons;
			electrons += change;
			neutrals -= change;
			// Don't allow neutrals dip below zero
			if (neutrals < 0.0) neutrals = 0.0;
			if (electrons > num_atoms) electrons = num_atoms;

			odeObj->rho[i + 1] = electrons;
		}

		odeObj->j_e[0] = j_e0;
		fv1 = j_e0;
		for (int i = 0; i < num_t - 1; i++)
		{
			fv2 = real(odeObj->ee_p[i + 1] + odeObj->ee_m[i + 1]);
			current_change = ht * (pow(charge_e, 2) / mass_e) * electrons * (fv1 + fv2) * 0.5 - ht * ve * current;
			fv1 = fv2;
			current += current_change;
			odeObj->j_e[i + 1] = current;
		}

		//if (z == distanceSourceToSample)
		

		// I THINK.. this takes the Current j_e and forewardFFTs it into the array called jhat which is then used in the Integrate()
		//applyWindow(odeObj->j_e);
		fftw_execute(odeObj->j_ffft);
		normalizeFFT(odeObj->jhat, fftnorm);
	}
	else if (odeObj->doPlasmaCalc == 2) {

		double ht = domain_t / num_td;
		double neutrals = num_atoms - rho_0;                          // Neutral particles
		double electrons = rho_0;                      // background Electrons
		double change = 0.0;    
		double current = 0.0;                                       // Current to be exported to UPPE
		double current_change = 0.0e0;                              // Current change for differential equation
		double ve = 0.0; //1/tauCollision;
		double fv1 = 0.0e0, fv2 = 0.0e0;

		odeObj->rho[0] = electrons;
		for (int i = 0; i < num_t - 1; i++)
		{
			double fieldIntensity = ( pow(real(odeObj->ee_p[i] + odeObj->ee_m[i]),2) ) / Znaught ;  // MIRO real+real
			change = neutrals * odeObj->mpi_sigmaK * ht * pow(fieldIntensity, odeObj->mpi_k);
			electrons += change;
			neutrals -= change;
			// Don't allow neutrals dip below zero
			if (neutrals < 0.0) neutrals = 0.0;
			if (electrons > num_atoms) electrons = num_atoms;

			odeObj->rho[i + 1] = electrons;
		}

		odeObj->j_e[0] = j_e0;
		fv1 = j_e0;
		for (int i = 0; i < num_t - 1; i++)
		{
			//odeObj->j_e[i + 1] = (1.0 - ht / tauCollision)*odeObj->j_e[i] + ht * pow(charge_e, 2) / mass_e * odeObj->rho[i] * real(odeObj->eFieldPlus[i] + odeObj->eFieldMinus[i]);
			// Andrew original
			//odeObj->j_e[i + 1] = odeObj->j_e[i] * exp(-1.0 * ht / tauCollision) + pow(charge_e, 2) / mass_e * 0.5 * ht * (odeObj->rho[i] * real(odeObj->ee_p[i] + odeObj->ee_m[i]) * exp(-1.0 * ht / tauCollision) + odeObj->rho[i + 1] * real(odeObj->ee_p[i + 1] + odeObj->ee_m[i + 1]));
			
			
			// CALCULATING THE CURRENT for UPPE (based on Ewan's notes Jan-24 2018)
			//  First order method	
			//	current_change = delta_t*(pow(cnst_e,2)/cnst_me)*electrons*real(aptr[t])-delta_t*ve*current;	
			// Second order method (Kolja)
			fv2 = real(odeObj->ee_p[i + 1] + odeObj->ee_m[i + 1]);
			current_change = ht * (pow(charge_e, 2) / mass_e) * electrons * (fv1 + fv2) * 0.5 - ht * ve * current;
			fv1 = fv2;
			current += current_change;
			odeObj->j_e[i + 1] = current;
		}

		//if (z == distanceSourceToSample)
		

		// I THINK.. this takes the Current j_e and forewardFFTs it into the array called jhat which is then used in the Integrate()
		//applyWindow(odeObj->j_e);
		fftw_execute(odeObj->j_ffft);
		normalizeFFT(odeObj->jhat, fftnorm);

	}
	else if (odeObj->doPlasmaCalc == 1){
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
		double potentialFrac = odeObj->ionE / U_H;
		double potFrac52 = pow(potentialFrac, 5.0/2.0);
		double potFrac32 = pow(potentialFrac, 3.0/2.0);
		double wQST, eField, eFieldRatio;

		odeObj->rho[0] = electrons;
		for (int i = 0; i < num_t - 1; i++)
		{
			eField = real(odeObj->ee_p[i] + odeObj->ee_m[i]);
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

			odeObj->rho[i + 1] = electrons;


			/* if (eField != 0) {
				odeObj->j_e[i + 1] = odeObj->ionE * wQST * neutrals / eField;
			}
			else {
				odeObj->j_e[i + 1] = 0.0;
			} */
		}

		odeObj->j_e[0] = j_e0;
		fv1 = j_e0;
		for (int i = 0; i < num_t - 1; i++)
		{
			//odeObj->j_e[i + 1] = (1.0 - ht / tauCollision)*odeObj->j_e[i] + ht * pow(charge_e, 2) / mass_e * odeObj->rho[i] * real(odeObj->eFieldPlus[i] + odeObj->eFieldMinus[i]);
			// Andrew original
			//odeObj->j_e[i + 1] = odeObj->j_e[i] * exp(-1.0 * ht / tauCollision) + pow(charge_e, 2) / mass_e * 0.5 * ht * (odeObj->rho[i] * real(odeObj->ee_p[i] + odeObj->ee_m[i]) * exp(-1.0 * ht / tauCollision) + odeObj->rho[i + 1] * real(odeObj->ee_p[i + 1] + odeObj->ee_m[i + 1]));
			
			
			// CALCULATING THE CURRENT for UPPE (based on Ewan's notes Jan-24 2018)
			//  First order method	
			//	current_change = delta_t*(pow(cnst_e,2)/cnst_me)*electrons*real(aptr[t])-delta_t*ve*current;	
			// Second order method (Kolja)
			fv2 = real(odeObj->ee_p[i+1] + odeObj->ee_m[i+1]);
			current_change = ht * (pow(charge_e, 2) / mass_e) * electrons * (fv1 + fv2) * 0.5 - ht * ve * real(odeObj->j_e[i]);
			fv1 = fv2;
			odeObj->j_e[i + 1] = odeObj->j_e[i] + current_change;
		}
		

		//applyWindow(odeObj->j_e);
		fftw_execute(odeObj->j_ffft);
		normalizeFFT(odeObj->jhat, fftnorm);
	}
	else if (odeObj->doPlasmaCalc == 0){
		#pragma omp parallel for
		for (int i = 0; i < num_t; i++)
		{
			odeObj->jhat[i] = 0.0;
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
		//complex<double> deltazA = (1.0i*pow(odeObj->omega[i], 2) / (2.0*(odeObj->k[i])*clightSquared)*odeObj->p_nl[i] + odeObj->omega[i] / (2.0*(odeObj->k[i])*clightSquared*epsilon_0)*odeObj->jhat[i])*exp(-1.0i*real(odeObj->k[i])*z)*exp(-1.0*abs(imag(odeObj->k[i]))*z);
		complex<double> deltazA = (1.0i*pow(odeObj->omega[i], 2) / (2.0*(odeObj->k[i])*clightSquared)*odeObj->p_nl[i] + odeObj->omega[i] / (2.0*(odeObj->k[i])*clightSquared*epsilon_0)*odeObj->jhat[i])*exp(-1.0i*(odeObj->k[i])*z);
		f[i] = real(deltazA);
		f[i + num_tOver2 + 1] = imag(deltazA);
		//deltazA = -(1.0i*pow(odeObj->omega[i], 2) / (2.0*(odeObj->k[i])*clightSquared)*odeObj->p_nl[i] + odeObj->omega[i] / (2.0*(odeObj->k[i])*clightSquared*epsilon_0)*odeObj->jhat[i])*exp(1.0i*real(odeObj->k[i])*z)*exp(-1.0*abs(imag(odeObj->k[i]))*z);
		deltazA = -(1.0i*pow(odeObj->omega[i], 2) / (2.0*(odeObj->k[i])*clightSquared)*odeObj->p_nl[i] + odeObj->omega[i] / (2.0*(odeObj->k[i])*clightSquared*epsilon_0)*odeObj->jhat[i])*exp(1.0i*(odeObj->k[i])*z);
		f[i + num_t + 2] = real(deltazA);
		f[i + 3 * num_tOver2 + 3] = imag(deltazA); 

	}
	

	return GSL_SUCCESS;
}



void integrate(double z, double zStep, ODEParams *odeObj, double*y, complex<double>*integral) {
	for (int i = 0; i < numActiveOmega; i++)
	{
		if (i < freqLowerCutoff || i > freqUpperCutoff)
		{
			integral[i] = 0.0;
		}
		else {
			integral[i] += -(1.0i*pow(odeObj->omega[i], 2) / (2.0*(odeObj->k[i])*pow(cLight, 2))*odeObj->p_nl[i] + odeObj->omega[i] / (2.0*(odeObj->k[i])*pow(cLight, 2)*epsilon_0)*odeObj->jhat[i])* exp(-1.0i*real(odeObj->k[i]) * z)*exp(-1.0*abs(imag(odeObj->k[i]))*z)*zStep;
		}
	}
	return;
}

