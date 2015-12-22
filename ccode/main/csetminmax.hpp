#ifndef CSETMINMAX_H
#define CSETMINMAX_H

#include <vector>

std::vector<double> adjustToLimits(std::vector<double>, double, double, double);
double adjustToLimits(double, double, double, double);
std::vector<double> cSetMinMax(std::vector<double>, double, double);
double cSetMinMax(double, double, double);

#endif