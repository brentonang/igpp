#ifndef CPATHTOPOSRTM1_H
#define CPATHTOPOSRTM1_H

#include <vector>
#include <tuple>

const double db_res = 2.7778e-4; // 1/3600
const double hdb_res = db_res/2;
const int SRTM_SIZE = 1201;
const int row = 500;
const int col = 1000;

std::vector<int> returnMin(std::vector<int>, int);
// std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> Cpath_topoSRTM1(std::vector<double>, std::vector<double>);
std::pair<std::vector<double>, std::vector<double>> Cpath_topoSRTM1(std::vector<double>, std::vector<double>);
#endif