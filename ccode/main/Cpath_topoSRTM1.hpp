#ifndef CPATHTOPOSRTM1_H
#define CPATHTOPOSRTM1_H

#include <vector>
#include <tuple>

const double pi = 3.14159;
const double rad2deg = 180/pi;
const double deg2rad = pi/180;

std::vector<int> returnMin(std::vector<int>, int);
std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> Cpath_topoSRTM1(std::vector<double>, std::vector<double>);
#endif