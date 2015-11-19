#include <iostream>

using namespace std;

const float srclat = 41.131, srclon = 360 - 112.8965;	// source location
const int rinc = 100;									// range increment (m)
const float azimuth = 70;								// azimuth cw from N (degrees)
const int range = 160000;								// maximum range (m)
const int maxrange = 130000;

/*************************************************************************************
 * 1) Get vector of ranges along the path (maxrange <= 130km)
*************************************************************************************/
int max_range = min(maxrange, range);
int ranges = 
int main () {
	cout << "Hello World " << srclat << " " << srclon << endl;
	return 0;
}