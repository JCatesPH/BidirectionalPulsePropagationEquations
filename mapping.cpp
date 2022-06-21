#include "Simulation.h"

int dAdz(double z, const double y[], double f[], void *params);
void integrate(double z, double zStep, ODEParams *odeObj, double*y, complex<double>*integral);
void updateGuess(double *ynew, complex<double> *sLeft, const gsl_vector *guessAm, RootParams *rootObj);

double min_mapG(const gsl_vector *ym_guess, void *simPtr) {
    //rootparam_type *rparams = reinterpret_cast<rootparam_type*>(rootparams);
    Simulation* simObj = reinterpret_cast<Simulation*>(simPtr);
    RootParams *rootObj = reinterpret_cast<RootParams*>(simObj->getRootParams());

    complex<double>* sLeft = simObj->getSourceLHS();
    complex<double>* sRight = simObj->getSourceRHS();

    int numOmeg = simObj->getNumOmega();
    int lowerFreq = simObj->getLowerFreqIndex();
    vector<double> pointMons = simObj->getZmonLocs();

	// Initialize GSL ODE objects
	int arrSize, GSLerrorFlag, foundNaN = 0;
	if (simObj->getNumDimsMinusOne() == 1) {
		arrSize = simObj->getNumOmX();
	}
	else {
		arrSize = numOmeg;
	}
	//const gsl_odeiv2_step_type * stepType = gsl_odeiv2_step_rkf45;
	const gsl_odeiv2_step_type * stepType = gsl_odeiv2_step_rk8pd;
	gsl_odeiv2_step *gslStep = gsl_odeiv2_step_alloc(stepType, 4 * arrSize);
	gsl_odeiv2_evolve *gslEvolve = gsl_odeiv2_evolve_alloc(4 * arrSize);
	//gsl_odeiv2_control * gslControl = gsl_odeiv2_control_standard_new(simObj->getODE_epsAbs(), simObj->getODE_epsRel(), 0.9, 0.1); // 3rd and 4th parameter set scaling factors of y(t) and y'(t) respectively
	gsl_odeiv2_control * gslControl = gsl_odeiv2_control_y_new(simObj->getODE_epsAbs(), simObj->getODE_epsRel());
    gsl_odeiv2_system sys = { dAdz, NULL, (size_t)(4 * arrSize), rootObj->getODEparams()};
	
	gsl_odeiv2_evolve_reset(gslEvolve);

	double *yloc = (rootObj->getODEparams())->y;

    double nonlinear_time_initial, nonlinear_time;

	int numzReports, numZsteps = 0;
	double zRight, zStepSize;
	double zPosition = 0.0;
	vector<double> endPoints;

    nonlinear_time_initial = omp_get_wtime();

	// Reset the left source and update guess
	updateGuess(yloc, sLeft, ym_guess, rootObj);
	
	if (rootObj->getOutParam() == 1) {
		write_out_eFieldAndSpectrumAtZlocation(rootObj->getItNum(), 
            0, yloc, 0.0, rootObj->getODEparams()->ee_m, 
            simObj->getStructure()->m_layers.front().getMaterial().getK(), 
            rootObj->getODEparams()->em_b);
	}
    //cout << "  Going FORWARD through layers" << endl;
    for (std::list<Layer>::iterator lit = simObj->getStructure()->m_layers.begin(); lit != simObj->getStructure()->m_layers.end(); ++lit) {
        // Skip the LHS layer and the RHS layers
        if (lit->getLowSideBoundary() != NULL && lit->getHiSideBoundary() != NULL)
        {
            boundary(lit->getLowSideBoundary()->m_zPos, 
				lit->getLowSideBoundary()->lowSideLayer->getMaterial().getK(), 
				lit->getLowSideBoundary()->hiSideLayer->getMaterial().getK(), 
				yloc);
		}
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

            for (std::size_t atZpos = 0; atZpos < pointMons.size(); ++atZpos) {
                if (pointMons[atZpos] > lit->getStartZpos() && pointMons[atZpos] < lit->getEndZpos())
                {
                    endPoints.push_back(pointMons[atZpos]);
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

					for (int k = 0; k < arrSize; k++){
						if(yloc[k] != yloc[k]) foundNaN++;
						if(yloc[k + arrSize] != yloc[k + arrSize]) foundNaN++;
						if(yloc[k + 2 * arrSize] != yloc[k + 2 * arrSize]) foundNaN++;
						if(yloc[k + 3 * arrSize] != yloc[k + 3 * arrSize]) foundNaN++;
					}

					if (foundNaN >= 1) {
						printf("WARNING: found NaN\n");
						printf("   Info: While taking ode step at z = %.5e\n", zPosition);
						raise(SIGABRT);
					}
                    
					if (GSLerrorFlag == GSL_SUCCESS) {
                        numZsteps++;

						if (rootObj->getIntCond() != 0) {
							integrate(zPosition, zStepSize, rootObj->getODEparams(), yloc, rootObj->integral);
						}
                    }
                    else {
                        printf("error: ODE driver returned %d\n", GSLerrorFlag);
                        gsl_odeiv2_evolve_reset(gslEvolve);
                    }

                    nonlinear_time = omp_get_wtime() - nonlinear_time_initial;
                    if ((int)(nonlinear_time / 20) > numzReports && rootObj->getOutParam() == 1) {
                        printf("  I = %d, step = %d, z = %.8g, t = %d s\n", rootObj->getItNum(), numZsteps, zPosition, (int)nonlinear_time);
                        numzReports++;
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
            
        //}
        //  Finally do the lowside boundary of the last assumed Vacuum layer
        if (lit->getHiSideBoundary() == NULL) {
            boundary(lit->getLowSideBoundary()->m_zPos, 
				lit->getLowSideBoundary()->lowSideLayer->getMaterial().getK(), 
				lit->getLowSideBoundary()->hiSideLayer->getMaterial().getK(), 
				yloc);
        }
    }


	/* for (int k = 0; k < rootObj->getSizeRoot()/2; k++){
		gsl_vector_set(f, k, yloc[k + 2 * arrSize + freqLowerCutoff] - real(sourceRight[k + freqLowerCutoff]));
		gsl_vector_set(f, k + rootObj->getSizeRoot()/2, yloc[k + 3 * arrSize + freqLowerCutoff] - imag(sourceRight[k + freqLowerCutoff]));
	} */
	double f2norm = 0.0;
	for (int k = 0; k < rootObj->getSizeRoot()/2; k++){
		f2norm += pow(yloc[k + 2 * numOmeg + lowerFreq] - real(sRight[k + lowerFreq]), 2);
		f2norm += pow(yloc[k + 3 * numOmeg + lowerFreq] - imag(sRight[k + lowerFreq]), 2);
	}
	f2norm = sqrt(f2norm);
	rootObj->setGnorm(f2norm);


	if (rootObj->getOutParam() == 1) {
		write_out_eFieldAndSpectrumAtZlocation(rootObj->getItNum(), 
        1, yloc, simObj->getStructure()->getThickness(), rootObj->getODEparams()->ee_p, 
        simObj->getStructure()->m_layers.back().getMaterial().getK(), 
        rootObj->getODEparams()->ep_b);
	}
    
	nonlinear_time = omp_get_wtime() - nonlinear_time_initial;
	//if (rootObj->output == 1) {
		//printf("Iteration %d completed in %.2f seconds with %d steps.\n", rootObj->itnum, nonlinear_time, numZsteps);
		//cout << "Iteration " << rootObj->itnum <<  " completed in " <<  nonlinear_time << "seconds with" << numZsteps << "steps." << endl;
		//fflush(stdout);
	//}

	for (int k = 0; k < arrSize; k++){
		yloc[k + 2*arrSize] = real(sRight[k]);
		yloc[k + 3*arrSize] = imag(sRight[k]);
	}

	simObj->getStructure()->doBackwardPassThroughAllBoundaries(yloc);

	// Freeing ODE memory
	gsl_odeiv2_control_free(gslControl);
    gsl_odeiv2_evolve_free(gslEvolve);
    gsl_odeiv2_step_free(gslStep);

	for (int k = 0; k < arrSize; k++){
		if(yloc[k] != yloc[k]) foundNaN++;
		if(yloc[k + arrSize] != yloc[k + arrSize]) foundNaN++;
		if(yloc[k + 2 * arrSize] != yloc[k + 2 * arrSize]) foundNaN++;
		if(yloc[k + 3 * arrSize] != yloc[k + 3 * arrSize]) foundNaN++;
	}

	/* if (foundNaN >= 1) 
		return GSL_EBADFUNC;
	else
    	return GSL_SUCCESS; */
	return f2norm;
}
