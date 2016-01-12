#include "Cpath_topoSRTM1.hpp"
#include "csetminmax.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

vector<int> returnMin(vector<int> firstVector, int secondValue) {
    vector<int> returnVector;
    for (int i = 0; i < firstVector.size(); i++) returnVector[i] = (firstVector[i] < secondValue) ? firstVector[i] : secondValue;
    return returnVector;
}

tuple<vector<double>, vector<double>, vector<double>> Cpath_topoSRTM1(vector<double> path_lat, vector<double> path_lon) {
    char* buffer;
    FILE *pFile;
    int i, status;
    long fileSize;
    double db_res, hdb_res, nbyters_per_lon;
    vector<double> glon(10801), glat(14401), olon(1302), olat(1302);
    vector<int> ilat(1302), ilon(1302), ilatplus1(1302), ilonplus1(1302), minLat(1302), minLon(1302), offset(1302);
    vector<int> db_loc{39, 43, 247, 250};
    vector<int> db_size{14401, 10801};
    tuple<vector<double>, vector<double>, vector<double>> returnTuple;

    db_res = 1/3600;
    hdb_res = db_res/2;
    nbyters_per_lon = 2.0 * db_size[0];

    path_lon = cSetMinMax(path_lon, 0, 360);
    for (i = 0; i < path_lon.size(); i++) {
        ilat[i] = round((path_lat[i] - db_loc[0])/db_res);
        ilon[i] = round((path_lon[i] - db_loc[2])/db_res);
        offset[i] = ilon[i] * nbyters_per_lon + ilat[i] * 2;
    }

    pFile = fopen("N39W113_4X3.hgt", "rb");
    if(pFile == NULL) {
        cout << "Error opening file" << endl;
        vector<double> zeros1(path_lat.size());
        vector<double> zeros2(path_lat.size());
        vector<double> zeros3(path_lat.size());
        tuple<vector<double>, vector<double>, vector<double>> zeroTuple;
        std::get<0>(zeroTuple) = zeros1;
        std::get<1>(zeroTuple) = zeros2;
        std::get<2>(zeroTuple) = zeros3;
        return zeroTuple;
    }
    vector<double> path_data(path_lon.size());
    fileSize = ftell(pFile);
    rewind(pFile);
    buffer = (char*) malloc (sizeof(char)*fileSize);
    for (i = 0; i < path_lon.size(); i++) {
        status = fseek(pFile, offset[i], SEEK_SET); // SEEK_SET == bof
        path_data[i] = fread(buffer, sizeof(int) * 2, 1, pFile); // fread(ptr, size, count, stream)
    }

    fclose(pFile);
    free(buffer);
    for (i = 0; i < db_size[1]; i++) {
        glon[i] = (i + 1) * db_res - hdb_res;
        for (int j = 0; j < db_size[0]; i++) glat[j] = 90 - (j + 1) * db_res + hdb_res;
        for (int k = 0; k < ilat.size(); k++) ilatplus1[k] = ilat[k] + 1;
        for (int l = 0; l < ilon.size(); l++) ilonplus1[l] = ilon[l] + 1;
        minLon = returnMin(ilonplus1, db_size[1]);
        minLat = returnMin(ilatplus1, db_size[0]);
        for (int m = 0; m < minLon.size(); m++) olon[m] = glon[minLon[m]] - hdb_res;
        for (int n = 0; n < minLat.size(); n++) olat[n] = glat[minLat[n]] + hdb_res;
    }
    std::get<0>(returnTuple) = path_data;
    std::get<1>(returnTuple) = olat;
    std::get<2>(returnTuple) = olon;
    return returnTuple;
}

// int main() {
//     srand(time(0));
//     int randomVariable1 = rand();
//     int randomVariable2 = rand();
//     vector<double> vector1(1302, randomVariable1);
//     vector<double> vector2(1302, randomVariable2);
//     tuple<vector<double>, vector<double>, vector<double>> returnTuple;
//     returnTuple = Cpath_topoSRTM1(vector1, vector2);
//     vector<double> firstValue = std::get<0>(returnTuple);
//     vector<double> secondValue = std::get<1>(returnTuple);
//     vector<double> thirdValue = std::get<2>(returnTuple);
//     cout << "The value of the first tuple is: ";
//     for (vector<double>::iterator it = firstValue.begin(); it != firstValue.end(); it++) cout << *it << " ";
//     cout << endl;
//     cout << "The value of the second tuple is: ";
//     for (vector<double>::iterator it = secondValue.begin(); it != secondValue.end(); it++) cout << *it << " ";    
//     cout << endl;
//     cout << "The value of the third tuple is: ";
//     for (vector<double>::iterator it = thirdValue.begin(); it != thirdValue.end(); it++) cout << *it << " ";    
// }