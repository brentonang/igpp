#include <iostream>
#include <stdlib.h>
#include <vector>

using namespace std;

int main(int argc, char* argv[]) {
    int i, j;
    float srclat, srclon, rinc, azimuth, range, maxrange, max_range;
    vector<float> flatlon, flatlat, plon, ranges;
    if (argc != 7) {
        // cout << "Invalid number of command line arguments" << endl;
    } else {
        srclat = atof(argv[1]);
        srclon = atof(argv[2]);
        rinc = atof(argv[3]);
        azimuth = atof(argv[4]);
        range = atof(argv[5]);
        maxrange = atof(argv[6]);
    }
    max_range = min(maxrange, range);

    i= 0;
    j = 0;
    while (j <= max_range+rinc) {
        ranges[i] += j;
        j += 100;
    } 
    
    if (azimuth < 0 || azimuth > 180) {
        cout << "ERROR: Azimuth must be eastward" << endl;
        return 0;
    }

    vector<float> azm(ranges.size(), 1*azimuth);

    // cget_latlon(srclat, srclon, ranges/1852, azm);
    // plon = csetminmax(plon, 0, 360);

    return 0;
}