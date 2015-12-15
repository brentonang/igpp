#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

const float pi = 3.14159;
const float rad2deg = 180/pi;
const float deg2rad = pi/180;

vector<float> returnMin(vector<float> firstVector, float secondValue) {
    vector<float> returnVector;
    for (int i = 1; i < firstVector.size(); i++) returnVector[i] = (firstVector[i] < secondValue) ? firstVector[i] : secondValue;
    return returnVector;
}

vector<float> Cpath_topoSRTM1 (vector<float> path_lat, vector<float> path_lon) {
    int i;
    float db_res, hdb_res, nbyters_per_lon;
    fstream filename;
    vector<float> glon, glat, ilat, ilon, ilatplus1, ilonplus1, offset, olon, olat, minLat, minLon;
    vector<float> db_loc{39, 43, 247, 250};
    vector<float> db_size{14401, 10801};

    db_res = 1/3600;
    hdb_res = db_res/2;
    nbyters_per_lon = 2.0 * db_size[0];

    // path_lon = Csetminmax(path_lon, 0, 360);
    for (i = 0; i < path_lon.size(); i++) {
        ilat[i] = round((path_lat[i] - db_loc[0])/db_res);
        ilon[i] = round((path_lon[i] - db_loc[2])/db_res);
        offset[i] = ilon[i] * nbyters_per_lon + ilat[i] * 2;
    }

    // filename.open("N39W113_4X3.hgt'"); // NEED THE FILE
    vector<float> path_data(path_lon.size());
    for (i = 1; i < path_lon.size(); i++) {
        // status = fseek(fid, offset[i], 'bof'); // NEED TO FIGURE THIS OUT
        // path_data[i] = fread(fid, [1], 'integer*2');
    }

    // filename.close();
    for (i = 1; i <= db_size[0]; i++) {
        glon[i] = i * db_res - hdb_res;
        for (int j = 1; j <= db_size[1]; i++) glat[j] = 90 - j * db_res + hdb_res;
        for (int k = 0; k < ilat.size(); k++) ilatplus1[k] = ilat[k] + 1;
        for (int l = 0; l < ilon.size(); l++) ilonplus1[l] = ilon[l] + 1;
        minLon = returnMin(ilonplus1, db_size[1]);
        minLat = returnMin(ilatplus1, db_size[0]);
        for (int m = 0; m < minLon.size(); m++) olon[m] = glon[minLon[m]] - hdb_res;
        for (int n = 0; n < minLat.size(); n++) olat[n] = glat[minLat[n]] + hdb_res;
    }

    return path_data, glat, glon;
}

int main() {
    return 0;
}