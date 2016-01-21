#ifndef CGETLATLON_H
#define CGETLATLON_H

#include <cmath>
#include <vector>

const double pi = 3.14159;
const double deg2rad = pi/180;
const double rad2deg = 180/pi;
const double halfpi = pi/2;
// const double dtol = 1*pow(10, -12);
const double dtol = 1.0e-12;
const double meritol = 1*pow(10, -10);

std::vector<double> elementWiseMultiplication(std::vector<double>, double);
std::vector<double> elementWiseDivision(std::vector<double>, double);
std::pair<std::vector<double>, std::vector<double>> cGetLatLon(double, double, std::vector<double>, std::vector<double>, double, double);

#endif