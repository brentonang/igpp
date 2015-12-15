/**************************
 * RESEARCH OPEN SOURCE CODE PERMISSIONS
 *
 *************************/


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

// blat and blon are single valued variables not vector or array
vector<float> cGetLatLon (int argc, vector<float> blat, vector<float> blon, vector<float> ranges, vector<float> azm, float major_axis, float esquared) {
	float ellipse, onef, f, f4;
	int i, j, len = ranges.size();
	vector<float> c1, c2, S, phi1, phi2, lam1, lam2, al12, al21, signS, sina12, cosa12, th1, costh1, sinth1, merid, D, M, N, P, s1, d, u, V;
	vector<float> sind, X, ds, ss, cosds, sinds, de, denom, plon, plat;
	// if (argc < 5) {
	// 	major_axis =  6378206.4;
	// 	esquared = 0.006768658;  
	// } OPTION TO PUT IN ELLIPSOID

	if (major_axis == 0) {
		major_axis = 6378206.4;
		esquared = 0;
	} // OPTION TO PUT IN CIRCLE

	// do {
	//  	cout << "Press the Enter Key to continue";
	// } while (std::cin.get() != '\n');
    // KEYBOARD equivalent of C++ where you can check variable values

	if (azm.size() != len) {
		cout << "Ranges and AZM must be the same length" << endl;
		vector<float> zeros(len);
		return zeros;
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

	// al12 = cadjlon(al12); RUN FIRST IN MATLAB AND SEE WHAT VALUE IS AND CHECK AGAINST CSETMINMAX
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
		if (abs(sina12[i]) < meritol) {
			merid[i] = 1;
			if (ellipse == 1) {
				c1[i] = 0;
				c2[i] = f4;
				D[i] = (1 - c2[i]);
				D[i] *= D[i];
				P[i] = c2[i]/D[i];
			}
		}
		else {
			merid[i] = 0;
			if (ellipse == 1) {
				c1[i] = f * M[i];
				c2[i] = f4 * ((1 - M[i]) * M[i]);
				D[i] = (1 - c2[i]) * (1 - c2[i] - c1[i] * M[i]);
				P[i] = ((1 + 5 * c1[i]) * M[i]) * (c2[i]/D[i]);
			}
		}
	}

	for (i = 1; i < len; i++) {
		if (merid[i] == 1) {
			sina12[i] = 0;
			for (j = 1; j < th1.size(); j++) s1[i] = halfpi - th1[j]; // DOUBLE CHECK LOGIC
			if (abs(al12[i]) < halfpi) cosa12[i] = 1;
			else cosa12[i] = -1;
			M[i] = 0;
		} else {
			s1[i] = (abs(M[i]) >= 1) ? 0 : acos(M[i]);
			s1[i] = sinth1[i]/sin(s1[i]);
			s1[i] = (abs(s1[i]) >= 1) ? 0 : acos(s1[i]);
			for (i = 0; i < sina12.size(); i++) {
				M[i] = costh1[i] * sina12[i];
			}
		}
	}

	transform(costh1.begin(), costh1.end(), cosa12.begin(), N.begin(), std::multiplies<float>()); // TEST CORRECTNESS

	if (ellipse == 1) {
		vector<float> tempVector;
		tempVector = elementWiseMult(D, major_axis);
		transform(S.begin(), S.end(), tempVector.begin(), d.begin(), std::divides<float>());
		for (i = 0; i < al12.size(); i++ ) {
			if (abs(al12[i]) > halfpi) d[i] = -d[i];
		}
		for (j = 0; j < s1.size(); j++) { // right upper limit for j?
			u[j] = 2 * (s1[j] - d[j]);
			V[j] = cos(u[j] + d[j]);
			sind[j] = sin(d[j]);
			X[j] = c2[j] * c2[j] * sind[j] * cos(d[j]) * (2 * V[j] * V[j] - 1);
			ds[j] = d[j] + X[j] - 2 * P[j] * V[j] * (1 - 2 * P[j] * cos(u[j])) * sind[j];
			ss[j] = s1[j] + s1[j] - ds[j];
		}
	} else {
		for (j = 0; j < s1.size(); j++) {
			ds[j] = S[j] / major_axis;
		}
		for (i = 0; i < al12.size(); i++) {
			if (abs(al12[i]) > halfpi) ds[i] = -ds[i];
		}
	}

	for (i = 0; i < ds.size(); i++) {
		cosds[i] = cos(ds[i]);
		sinds[i] = sin(ds[i]);
	}
	for (i = 0; i < al12.size(); i++) {
		if (abs(al12[i]) > halfpi) sinds[i] = -sinds[i];
	}
	for (i = 0; i < sinth1.size(); i++) al21[i] = N[i] * cosds[i] - sinth1[i] * sinds[i];

	for (j = 1; j < len; j++) {
		if (merid[j] == 1) {
			phi2[j] = atan(tan(halfpi + s1[j] - ds[j])/onef);
			if (al21[j] > 0) {
				al21[j] = pi;
				if (signS[j] == 1) de[j] = pi;
				else {
					phi2[j] = -phi2[j];
					de[j] = 0;
				}
			} else {
				al21[j] = 0;
				if (signS[j] == 1) {
					phi2[j] = -phi2[j];
					de[j] = 0;
				} else de[j] = pi;
			}
		} else {
			al21[j] = atan(M[j]/al21[j]);
			if (al21[j] > 0) al21[j] += pi;
			if (al21[j] < 0) al21[j] -= pi;
			// al21[j] = cadjlon(al21[j]); // -pi to pi
			denom[j] = (ellipse == 1) ? onef * M[j] : M[j];
			phi2[j] = atan(-(sinth1[j] * cosds[j] + N[j] * sinds[j]) * sin(al21[j])/denom[j]);
			de[j] = atan2(sinds[j] * sina12[j], (costh1[j] * cosds[j] - sinth1[j] * sinds[j] * cosa12[j]);
			if (ellipse == 1) {
				if (signS[j] == 1) de[j] = de[j] + c1[j] * ((1 - c2[j]) * ds[j] + c2[j] * sinds[j] * cos(ss[j]));
				else de[j] = de[j] - c1[j] * ((1 - c2[j]) * ds[j] - c2[j] * sinds[j] * cos(ss[j]));
			}
		}
	}
	// lam2 = cadjlon(lam1 + de);
	lam1 = elementWiseMult(lam1, rad2deg);
	phi1 = elementWiseMult(phi1, rad2deg);
	al12 = elementWiseMult(al12, rad2deg);
	plon = elementWiseMult(lam2, rad2deg);
	plat = elementWiseMult(phi2, rad2deg);
	return plat, plon;
}

int main () {
	return 0;
}