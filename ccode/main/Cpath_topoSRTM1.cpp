#include "Cpath_topoSRTM1.hpp"
#include "csetminmax.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <tuple>
#include <vector>

using namespace std;

vector<double> returnMin(vector<double> firstVector, double secondValue) {
    vector<double> returnVector;
    for (int i = 1; i < firstVector.size(); i++) returnVector[i] = (firstVector[i] < secondValue) ? firstVector[i] : secondValue;
    return returnVector;
}

tuple<vector<double>, vector<double>, vector<double>> Cpath_topoSRTM1(vector<double> path_lat, vector<double> path_lon) {
    int i;
    double db_res, hdb_res, nbyters_per_lon;
    vector<double> glon, glat, ilat, ilon, ilatplus1, ilonplus1, offset, olon, olat, minLat, minLon;
    vector<double> db_loc{39, 43, 247, 250};
    vector<double> db_size{14401, 10801};
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

    ifstream file("N39W113_4X3.hgt", std::ios::in|std::ios::binary);
    if(!file) {
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
    for (i = 0; i < path_lon.size(); i++) {
        // status = fseek(fid, offset[i], 'bof'); // NEED TO FIGURE THIS OUT
        // path_data[i] = fread(fid, [1], 'integer*2');
    }

    file.close();
    for (i = 1; i <= db_size[1]; i++) {
        glon[i] = i * db_res - hdb_res;
        for (int j = 1; j <= db_size[1]; i++) glat[j] = 90 - j * db_res + hdb_res;
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