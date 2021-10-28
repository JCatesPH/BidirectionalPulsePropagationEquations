#pragma once

using namespace std;

std::string get_current_dir();
char *readParmetersFileToBuffer(char *paramFilePath);
int getIntParameterValueByName(char *paramName);
int getStringParameterValueByName(char *paramName, char *stringBuffer);
double getDoubleParameterValueByName(char *paramName);
double vnorm(const double *v, const int n);
