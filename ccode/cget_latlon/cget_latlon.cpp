#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

const float pi = 3.14159;
const float deg2rad = pi/180;
const float rad2deg = 180/pi;
const float halfpi = pi/2;
const float dtol = 1*pow(10, -12);
const float meritol = 1*pow(10, -10);

vector<float> elementWiseMult (vector<float> input, float value) {
	for (int i = 0; i < input.size(); i++) {
		input[i] = input[i] * value;
	}
	return input;
}

int cGetLatLon (int argc, vector<float> blat, vector<float> blon, vector<float> ranges, vector<float> azm, float major_axis, float esquared) {
	float ellipse, onef, f, f4, is1, is0, im1, im0;
	int i;
	vector<float> S, phi1, lam1, al12, signS, sina12, cosa12, th1, costh1, sinth1, merid, M, N;
	// if (argc < 5) {
	// 	major_axis =  6378206.4;
	// 	esquared = 0.006768658;  
	// }

	if (major_axis == 0) {
		major_axis = 6378206.4;
		esquared = 0;
	}

	do {
	 	cout << "Press the Enter Key to continue";
	} while (std::cin.get() != '\n');

	if (azm.size() != ranges.size()) {
		cout << "Ranges and AZM must be the same length" << endl;
		return 0;
	}

	S = elementWiseMult(ranges, 1852);
	phi1 = elementWiseMult(blat, deg2rad);
	lam1 = elementWiseMult(blon, deg2rad);
	al12 = elementWiseMult(azm, deg2rad);

	ellipse = ((esquared != 0.0) ? 1.0 : 0.0);

	onef = 1;
	f = 0;
	f4 = 0;
	if (ellipse == 1) {
		onef = sqrt(1 - esquared);
		f = 1 - onef;
		f4 = f/4;
	}

	// al12 = cadjlon(al12); SHOULD WE USE CSETMINMAX? UPLIM? LOWLIM?
	for (i = 0; i < al12.size(); i ++) { 	// DOUBLE CHECK
		sina12[i] = sin(al12[i]);
		cosa12[i] = cos(al12[i]);			// EARLY ASSIGNMENT OF SINAL12
		if (abs(al12[i]) > halfpi) signS[i] = 1;
		else signS[i] = 0;
	}

	th1 = phi1;
	if (ellipse == 1) {
		for (i = 0; i < phi1.size(); i++) {
			th1[i] = atan(onef*tan(phi1[i]));
		}
	}

	for (i = 0; i < th1.size(); i++) {
		costh1[i] = cos(th1[i]);
		sinth1[i] = sin(th1[i]);
	}

	for (i = 0; i < sina12.size(); i++) {
		if (abs(sina12[i]) < meritol) merid[i] = 1;
		else merid[i] = 0;
	}

	for (i = 1; i < ranges.size(); i++) {
		if (merid[i] == 1) {
			sina12[i] = 0;
			if (abs(al12[i]) < halfpi) cosa12[i] = 1;
			else cosa12[i] = 0;
			M[i] = 0;
		} else {
			for (i = 0; i < sina12.size(); i++) {
				M[i] = costh1[i] * sina12[i];
			}
		}
	}

	transform(costh1.begin(), costh1.end(), cosa12.begin(), N.begin(), std::multiplies<float>()); // TEST CORRECTNESS

	if (ellipse == 1){
		// LINE 115
	}
	return 0;
}

int main () {
	return 0;
}