#ifndef CGETLATLON_H
#define CGETLATLON_H

#include <vector>

class CGetLatLon {
public:
	std::vector<double> elementWiseMultiplication(std::vector<double>, double);
	std::vector<double> elementWiseDivision(std::vector<double>, double);
	std::vector<double> cGetLatLon(double, double, std::vector<double>, std::vector<double>, double, double);
	std::vector<double> cGetLatLon(double, double, std::vector<double>, std::vector<double>);
};

#endif