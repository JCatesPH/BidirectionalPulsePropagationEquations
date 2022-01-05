#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "BPPE.h"

using namespace std;

#ifdef _WIN64
    #include <direct.h>
    #define GetCurrentDir _getcwd
#else
    #include <unistd.h>
    #define GetCurrentDir getcwd
#endif


extern int VERBOSE;
extern char *paramFileBuffer;

std::string get_current_dir() {
	char buff[FILENAME_MAX]; //create string buffer to hold path
	GetCurrentDir(buff, FILENAME_MAX);
	string current_working_dir(buff);
	return current_working_dir;
}

char *readParmetersFileToBuffer(char *paramFilePath) {
	char *buffer = NULL;
    char *line = NULL;

	long length=0;
	ssize_t nread;
	size_t len = 0;

	FILE *parmFile = fopen(paramFilePath, "r");
	if (VERBOSE > 3) { printf("Reading Input paramaterfile  %s \n", paramFilePath ); }	

	if (parmFile)
	{
	  fseek(parmFile, 0, SEEK_END);
	  length = ftell(parmFile);
	  fseek(parmFile, 0, SEEK_SET);
	  buffer = (char *)malloc(length+1);
	  if(buffer)
	  {
			// put an empty space in the buffer so it has a null charactter
			sprintf(buffer," ");
			while ((nread = getline(&line, &len, parmFile)) != -1) {
				// skip lines beginning wit # or SPACE or EMPTY
               if(!((strncmp(line, "#", 1) == 0) || (strncmp(line, " ", 1) == 0) || (strlen(line) == 1)))
               {
					if (VERBOSE > 8) {printf("Valid NAME=VALUE found: %s", line);}
					sprintf(buffer + strlen(buffer), "%s", line);			   
			   }           
           }
	  }
	  fclose (parmFile);
	  // put an null charactter at the end just in case
	  buffer[length]='\0';
	}
	else
	{
		printf("ERROR: couldnt open parameter FILE  %s \n",paramFilePath);
		exit(-1);
	}	
	return buffer;			
}


int getIntParameterValueByName(char *paramName) {
  char *p;
  int found=0, paramValue=0;
  // have to copy the buffer  strtok() function modifies the buffer it searches 
  char *buffer=(char *)malloc(sizeof(char)*strlen(paramFileBuffer));
  strcpy(buffer,paramFileBuffer);
	    
  const char delims[24] = ",:= \r\n\t";
  p = strtok(buffer,delims);
  while (p!= NULL && found==0)
  {
    if (VERBOSE > 8) { printf("TOKEN: found token %s \n", p);}
	if (strcmp(p,paramName) == 0)
	{
		found = 1;
	    p = strtok (NULL, delims);
		paramValue = atoi(p);
	    if (VERBOSE > 8) { printf("SUCESS: found %s its value is %d\n", paramName, paramValue);}
	}
    p = strtok (NULL, delims);
  }
  free(buffer);
  if (found == 0) { printf("ERROR: Could not find a value %s in parameter file \n", paramName);exit(-1);}
  if (VERBOSE > 5) { printf("PARAMETER %s = %d\n", paramName, paramValue); }
  
  return paramValue;	
}

int getStringParameterValueByName(char *paramName, char *stringBuffer) {
  char *p;
  int found=0;  
  const char delims[24] = ",:= \r\n\t";
  
    // have to copy the buffer  strtok() function modifies the buffer it searches 
  char *buffer=(char *)malloc(sizeof(char)*strlen(paramFileBuffer));
  strcpy(buffer,paramFileBuffer);

  p = strtok(buffer,delims);
  while (p!= NULL && found==0)
  {
	if (strcmp(p,paramName) == 0)
	{
		found = 1;
	    p = strtok (NULL, delims);
	    strcpy(stringBuffer,p);
	    if (VERBOSE > 8) { printf("SUCESS: found %s its value is %s\n", paramName, p);}
	}
    p = strtok (NULL, delims);
  }
  free(buffer);

  if (found == 0) { printf("ERROR: Could not find a parameter    %s     in parameter file \n", paramName);exit(-1);}
  if (VERBOSE > 5) { printf("PARAMETER %s = %s\n", paramName, stringBuffer); }

  return found;	
}

double getDoubleParameterValueByName(char *paramName) {
  char *p;
  int found=0;
  double paramValue=0.0;
  // have to copy the buffer  strtok() function modifies the buffer it searches 
  char *buffer=(char *)malloc(sizeof(char)*strlen(paramFileBuffer));
  strcpy(buffer,paramFileBuffer);
	    
  const char delims[24] = ",:= \r\n\t";
  p = strtok(buffer,delims);
  while (p!= NULL && found==0)
  {
    if (VERBOSE > 8) { printf("TOKEN: found token %s \n", p);}
	if (strcmp(p,paramName) == 0)
	{
		found = 1;
	    p = strtok (NULL, delims);
		paramValue = atof(p);
	    if (VERBOSE > 8) { printf("SUCESS: found %s its value is %g\n", paramName, paramValue);}
	}
    p = strtok (NULL, delims);
  }
  free(buffer);
  if (found == 0) { printf("ERROR: Could not find a value %s in parameter file \n", paramName);exit(-1);}
  if (VERBOSE > 5) { printf("PARAMETER %s = %g\n", paramName, paramValue); }
  
  return paramValue;	
}


double vnorm(const double *v, const int n) {
	double norm = 0;
	for (int k = 0; k < n; k++) {
		norm += abs(v[k]);
	}
	return norm;
}


void fill_omg_k(double*omg, double*kx, MaterialDB &theMaterialDB) {

		for (int i = 0; i < num_t; i++)
		{
			if (i <= num_t / 2) {
				omg[i] = (2.0 * M_PI / domain_t)*i;
			}
			else {
				omg[i] = (2.0 * M_PI / domain_t)*((double)i - num_t);
			}
		}
		// COLM NEW this line replaces following comment block
		theMaterialDB.initAllMaterialKs(omg, numActiveOmega);

	return;
}

void initializeY(double *y, complex<double> *yp_init)
{
	for (int i = 0; i < numActiveOmega; i++)
	{
		y[i] = real(yp_init[i]);
		y[i + numActiveOmega] = imag(yp_init[i]);
		y[i + 2*numActiveOmega] = 0.0;
		y[i + 3*numActiveOmega] = 0.0;
	}

	for (int i = 0; i < freqLowerCutoff; i++) {
		y[i] = 0.0;
		y[i + numActiveOmega] = 0.0;
		y[i + 2*numActiveOmega] = 0.0;
		y[i + 3*numActiveOmega] = 0.0;
	}
	for (int i = freqUpperCutoff; i < numActiveOmega; i++) {
		y[i] = 0.0;
		y[i + numActiveOmega] = 0.0;
		y[i + 2*numActiveOmega] = 0.0;
		y[i + 3*numActiveOmega] = 0.0;
	}

}


void normalizeFFT(complex<double>* arr) {
	double c = sqrt((double)num_t);
	for (int i=0; i < num_t; i++) {
		arr[i] = arr[i] / c;
	}
}

double *window;
void createWindowFunc(double alpha){
	window = (double*)malloc(sizeof(double) * num_t);
	// Implements a Tukey window
	for (int i = 0; i < num_t/8; i++) {
		window[i] = 0.0;
	}
	for (int i=num_t/8; i < (int)((1+alpha)*num_t/8); i++) {
		window[i] = 0.5 * (1.0 - cos(4*M_PI*(i-num_t/4)/(alpha*num_t)));
	}
	for (int i=(int)((1+alpha)*num_t/8); i <= num_t/2; i++) {
		window[i] = 1.0;
	}
	for (int i=num_t/8; i <= num_t/2; i++) {
		window[num_t-i] = window[i];
	}
	for (int i = 7*num_t/8; i < num_t; i++) {
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

void applyWindow(complex<double>* arr){
	for (int i=0; i < num_t; i++) {
		arr[i] = arr[i] * window[i];
	}
} 