#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <string>
#include <list>
#include <omp.h>
#include <fftw3.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include "gsl/gsl_math.h"
#include "gsl/gsl_multimin.h"
#include "Utilities.h"
#include "Materials.h"
#include "Structure.h"
#include "gslParams.h"
#include "output.h"

// Includes for floating-point exception handling
#include <cfenv>
#include <csignal>

// Add functions for derivative and wrapper for optimization
void dGdA(const gsl_vector *Am_in, void *rootparams, gsl_vector *dG);
void GdG(const gsl_vector *Am_in, void *rootparams, double *G, gsl_vector *dG);

// Include struct for input pulse parameters
typedef struct {
	double A0;
	double relativeIntensity;
	double relativePhase;
	double pulseDuration;
	double omega0;
	double pulseWaist;
}pulseparam_type;

class Simulation
{
private:
    string m_simulationName;
    MaterialDB* m_materialDatabase = NULL;
    Structure *m_structure = NULL;
    ODEParams *m_ODEParams = NULL;
    RootParams *m_RootParams = NULL;

    int m_numDimensionsMinusOne = 0;
    int m_numT, m_numActiveOmega, m_numX, m_numOmX;
    double m_domT, m_domX;
    int m_freqLowerCutoff, m_freqUpperCutoff;
    double *m_omegaArray = NULL;
    complex<double>* m_sourceLeft = NULL;
    complex<double>* m_sourceRight = NULL;

    double m_ode_epsAbs, m_ode_epsRel;

    double m_rho0;

    vector<double> m_zMonLocs;

    double m_simTime_initial, m_simTime_current;
    int algType = -1; // Positive integer indicates algorithm type used
    int m_maxIter, m_outputInterval;
    double m_minInitStep, m_minStopCon, m_minTol;

    // Add functions for different types of algorithms
    void fdfMin_iterateBPPE(); // Integers 1-9
    void fMin_iterateBPPE();   // Integers 11-19
    void root_iterateBPPE();   // Integers 21-29

    // Add helper functions
    void generateTwoColorPulse(complex<double>* ee, fftw_plan e_f, complex<double>* source, pulseparam_type *pparams);
    void doLinearProblem(ODEParams *odeObj, complex<double> *sourceLeft, complex<double> *sourceRight, Structure& theStructure);
    void generateGuess(gsl_vector *u, RootParams *rootObj, ODEParams *odeObj);
    void updateGuess(double *ynew, complex<double> *sLeft, const gsl_vector *guessAm, RootParams *rootObj);
    

public:
    Simulation() { }
    Simulation(string theSimulationName) : m_simulationName(theSimulationName) { }

    // Add methods to set the material database and structure objects
    void setMaterialDatabase(MaterialDB* aMaterialDatabase) { m_materialDatabase = aMaterialDatabase; }
    MaterialDB* getMaterialDatabase() { return m_materialDatabase; }
    void setStructure(Structure* aStructure) { m_structure = aStructure; }
    Structure* getStructure() { return m_structure; }
    RootParams* getRootParams() { return m_RootParams; }
    int getNumDimsMinusOne() { return m_numDimensionsMinusOne; }
    int getNumOmega() { return m_numActiveOmega; }
    int getNumOmX() { return m_numOmX;}
    int getLowerFreqIndex() { return m_freqLowerCutoff; }
    int getUpperFreqIndex() { return m_freqUpperCutoff; }
    double getODE_epsAbs() { return m_ode_epsAbs; }
    double getODE_epsRel() { return m_ode_epsRel; }

    void readParamfile(char *inFile);

    void setZmonLocs(vector<double> pointMonLi) { m_zMonLocs = pointMonLi; }
    vector<double> getZmonLocs() { return m_zMonLocs; }

    // Methods for setting LHS and RHS sources
    void setSourceLHS(pulseparam_type *pparams) {
        complex<double>* eFieldPlus = (complex<double>*)malloc(sizeof(complex<double>)*m_numT); 
        fftw_plan eFieldPlusForwardFFT = fftw_plan_dft_1d(m_numT, reinterpret_cast<fftw_complex*>(&eFieldPlus[0]), reinterpret_cast<fftw_complex*>(&eFieldPlus[0]), FFTW_FORWARD, FFTW_WISDOM_TYPE);
        m_sourceLeft = (complex<double>*)malloc(sizeof(complex<double>) * m_numActiveOmega);
        generateTwoColorPulse(eFieldPlus, eFieldPlusForwardFFT, m_sourceLeft, pparams);
    }
    complex<double>* getSourceLHS() { return m_sourceLeft; }

    void setSourceRHS(pulseparam_type *pparams) {
        complex<double>* eFieldMinus = (complex<double>*)malloc(sizeof(complex<double>)*m_numT); 
        fftw_plan eFieldMinusForwardFFT = fftw_plan_dft_1d(m_numT, reinterpret_cast<fftw_complex*>(&eFieldMinus[0]), reinterpret_cast<fftw_complex*>(&eFieldMinus[0]), FFTW_FORWARD, FFTW_WISDOM_TYPE);
        m_sourceRight = (complex<double>*)malloc(sizeof(complex<double>) * m_numActiveOmega);
        generateTwoColorPulse(eFieldMinus, eFieldMinusForwardFFT, m_sourceRight, pparams);
    }
    complex<double>* getSourceRHS() { return m_sourceRight; }


};
