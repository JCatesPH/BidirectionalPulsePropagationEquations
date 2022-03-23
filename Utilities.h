#pragma once

using namespace std;


std::string get_current_dir();
char *readParmetersFileToBuffer(char *paramFilePath);
int getIntParameterValueByName(char *paramName);
int getStringParameterValueByName(char *paramName, char *stringBuffer);
double getDoubleParameterValueByName(char *paramName);
double vnorm(const double *v, const int n);
void fill_omg_k(double*omg, double*kx, MaterialDB &theMaterialDB);
//void initializeY(double *y, complex<double> *yp_init);
void normalizeFFT(complex<double>* arr, double c);

void createWindowFunc(double alpha);
void applyWindow(complex<double>* arr);
void floatingPointExceptions(int sig);
