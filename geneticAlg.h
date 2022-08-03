#pragma once
# include <cmath>
# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <iomanip>
# include <ctime>
# include <cstring>

#include "BPPE.h"

# define POPSIZE 100
# define MAXGENS 25
# define NVARS 1256
# define PXOVER 0.8
# define PMUTATION 0.2

struct genotype
{
  double gene[NVARS];
  double fitness;
  double upper[NVARS];
  double lower[NVARS];
  double rfitness;
  double cfitness;
};

void crossoverGA( int &seed );
void elitistGA( );
void evaluateG(RootParams *aRootParams);
int i4_uniform_ab ( int a, int b, int &seed );
void initializeGA(int &seed);
void keepBestGA( );
void mutateGA( int &seed );
double r8_uniform_ab ( double a, double b, int &seed );
void reportGA( int generation );
void selectorGA( int &seed );
void XoverGA( int one, int two, int &seed );
void getBestGA(gsl_vector *v);
