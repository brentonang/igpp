#include "csetminmax.hpp"
#include "cget_latlon.hpp"
#include "Cpath_topoSRTM1.hpp"
#include "CompareTopo_AlongPath.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <tuple>
#include <vector>
#include <utility>

using namespace std;

int main(int argc, char* argv[]) {
    int i = 0, j = 0;
    double azimuth, maxrange, max_range, range, rinc, scalar1, scalar2, srclat, srclon;
    vector<double> flatlat, flatlon, plat, plon, ranges, tempRanges, vector1, vector2;
    pair<vector<double>, vector<double>> tempValues;
    tuple<vector<double>, vector<double>, vector<double>> tempTuple;

    if(argc != 3) {
        cout << "Invalid number of command line arguments" << endl;
        return 0;
    } else {
        rinc = atof(argv[1]);
        if (rinc <= 30) {
            cout << "rinc should be greater than 30m" << endl;
            return -1;
        }
        azimuth = atof(argv[2]);
        if (azimuth < 0 || azimuth > 360) {
            cout << "azimuth should be between 0 and 360" << endl;
            return -1;
        }
        maxrange = atof(argv[3]);
        if (maxrange < 0 || maxrange > 150000) {
            cout << "maxrange should be between 0m and 150000m" << endl;
            return -1;
        }
    }

    max_range = min(maxrange, range);
    while(j <= max_range + rinc) {
        ranges[i] += j;
        j += rinc;
    } 
    
    if(azimuth < 0 || azimuth > 180) {
        cout << "ERROR: Azimuth must be eastward" << endl;
        return -1;
    }

    vector<double> azm(ranges.size(), 1 * azimuth);
    tempRanges = elementWiseDivision(ranges, 1852);
    tempValues = cGetLatLon(srclat, srclon, tempRanges, azm, 6378206.4, 0.006768658);
    plat = tempValues.first;
    plon = tempValues.second;
    plon = cSetMinMax(plon, 0, 360);

    scalar1 = srclon + sin(azimuth * deg2rad);
    for(i = 0; i < ranges.size(); i++) {
        vector1[i] = ranges[i] / (cos(srclat*deg2rad) * 111111);
    }
    flatlon = elementWiseMultiplication(vector1, scalar1);

    scalar2 = srclat + cos(azimuth * deg2rad);
    for(i = 0; i < ranges.size(); i++) {
        vector2[i] = ranges[i] / 111111;
    }
    flatlat = elementWiseMultiplication(vector2, scalar2);

    returnTuple = Cpath_topoSRTM1(plat, plon);

    

    return 0;
}