#include <algorithm>
#include <cmath>
#include <iostream>
#include "cget_latlon.hpp"
#include "csetminmax.hpp"

using namespace std;

const double pi = 3.14159;
const double deg2rad = pi/180;
const double rad2deg = 180/pi;
const double halfpi = pi/2;
const double dtol = 1*pow(10, -12);
const double meritol = 1*pow(10, -10);

vector<double> elementWiseMultiplication(vector<double> input, double value) {
	for (int i = 0; i < input.size(); i++) {
		input[i] = input[i] * value;
	}
	return input;
}

vector<double> elementWiseDivision(vector<double> input, double value) {
	for (int i = 0; i < input.size(); i++) {
		input[i] = input[i] / value;
	}
	return input;
}

vector<double> cGetLatLon(double blat, double blon, vector<double> ranges, vector<double> azm, double major_axis, double esquared) {
	double costh1, denom, ellipse, f, f4, lam1, onef, phi1, th1, sinth1;
	int i, j, len = ranges.size();
	vector<double> al12, al21, c1, c2, cosa12, cosds, D, d, de, ds, lam2, M, merid, N, P, phi2, plon, plat, S, s1, ss, signS, sina12, sind, sinds, u, V, X;

	if (major_axis == 0) {
		major_axis = 6378206.4;
		esquared = 0;
	} 

	if (azm.size() != len) {
		cout << "Ranges and AZM must be the same length" << endl;
		vector<double> zeros(len);
		return zeros;
	}

	al12 = elementWiseMultiplication(azm, deg2rad);
	S = elementWiseMultiplication(ranges, 1852);
	phi1 = blat * deg2rad;
	lam1 = blon * deg2rad;

	ellipse = (esquared != 0) ? 1 : 0;

	onef = 1;
	f = 0;
	f4 = 0;
	if (ellipse == 1) {
		onef = sqrt(1 - esquared);
		f = 1 - onef;
		f4 = f/4;
	}

	al12 = csetminmax(al12, 180, -180);
	for (i = 0; i < al12.size(); i ++) { 
		sina12[i] = sin(al12[i]);
		cosa12[i] = cos(al12[i]);			
		signS[i] = (abs(al12[i]) > halfpi) ? 1 : 0;
	}

	th1 = phi1;
	if (ellipse == 1) th1 = atan(onef * tan(phi1));
	costh1 = cos(th1);
	sinth1 = sin(th1);

	for (i = 0; i < sina12.size(); i++) {
		if (abs(sina12[i]) < meritol) {
			merid[i] = 1;
			if (ellipse == 1) {
				c1[i] = 0;
				c2[i] = f4;
				D[i] = 1 - c2[i];
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
				P[i] = (1 + 5 * c1[i] * M[i]) * (c2[i]/D[i]);
			}
		}
	}

	for (i = 1; i < len; i++) {
		if (merid[i] == 1) {
			sina12[i] = 0;
			cosa12[i] = (abs(al12[i]) < halfpi) ? 1 : -1;
			M[i] = 0;
			s1[i] = halfpi - th1;
		} else {
			M[i] = costh1 * sina12[i];
			s1[i] = (abs(M[i]) >= 1) ? 0 : acos(M[i]);
			s1[i] = sinth1/sin(s1[i]);
			s1[i] = (abs(s1[i]) >= 1) ? 0 : acos(s1[i]);
		}
	}

	N = elementWiseMultiplication(cosa12, costh1);
	if (ellipse == 1) {
		vector<double> tempVector;
		tempVector = elementWiseMultiplication(D, major_axis);
		transform(S.begin(), S.end(), tempVector.begin(), d.begin(), std::divides<double>());
		for (i = 0; i < al12.size(); i++ ) {
			if (abs(al12[i]) > halfpi) d[i] = -d[i];
		}
		for (j = 0; j < s1.size(); j++) {
			u[j] = 2 * (s1[j] - d[j]);
			V[j] = cos(u[j] + d[j]);
			sind[j] = sin(d[j]);
			X[j] = c2[j] * c2[j] * sind[j] * cos(d[j]) * (2 * V[j] * V[j] - 1);
			ds[j] = d[j] + X[j] - 2 * P[j] * V[j] * (1 - 2 * P[j] * cos(u[j])) * sind[j];
			ss[j] = s1[j] + s1[j] - ds[j];
		}
	} else {
		ds = elementWiseDivision(S, major_axis);
		// for (j = 0; j < s1.size(); j++) {
		// 	ds[j] = S[j] / major_axis;
		// }
		for (i = 0; i < al12.size(); i++) {
			if (abs(al12[i]) > halfpi) ds[i] = -ds[i];
		}
	}

	for (i = 0; i < ds.size(); i++) {
		cosds[i] = cos(ds[i]);
		sinds[i] = sin(ds[i]);
	}

	for (i = 0; i < al12.size(); i++) if (abs(al12[i]) > halfpi) sinds[i] = -sinds[i];
	for (i = 0; i < cosds.size(); i++) al21[i] = N[i] * cosds[i] - sinth1 * sinds[i];
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
			al21[j] = csetminmax(al21[j], 180, -180);
			denom = (ellipse == 1) ? onef * M[j] : M[j];
			phi2[j] = atan(-(sinth1 * cosds[j] + N[j] * sinds[j]) * sin(al21[j])/denom);
			de[j] = atan2(sinds[j] * sina12[j], costh1 * cosds[j] - sinth1 * sinds[j] * cosa12[j]);
			if (ellipse == 1) {
				de[j] = (signS[j] == 1) ? de[j] + c1[j] * ((1 - c2[j]) * ds[j] + c2[j] * sinds[j] * cos(ss[j])) : de[j] - c1[j] * ((1 - c2[j]) * ds[j] - c2[j] * sinds[j] * cos(ss[j]));
				// if (signS[j] == 1) de[j] = de[j] + c1[j] * ((1 - c2[j]) * ds[j] + c2[j] * sinds[j] * cos(ss[j]));
				// else de[j] = de[j] - c1[j] * ((1 - c2[j]) * ds[j] - c2[j] * sinds[j] * cos(ss[j]));
			}
		}
	}
	for(i = 0; i < de.size(); i++) {
		de[i] += lam1;
	}
	lam2 = csetminmax(de, 180, -180);
	lam1 *= rad2deg;
	phi1 *= rad2deg;
	al12 = elementWiseMultiplication(al12, rad2deg);
	plon = elementWiseMultiplication(lam2, rad2deg);
	plat = elementWiseMultiplication(phi2, rad2deg);
	return plat, plon;
}

vector<double> cGetLatLon(double blat, double blon, vector<double> ranges, vector<double> azm) {
	double major_axis = 6378206.4;
	double esquared = 0.006768658;
	cGetLatLon(blat, blon, ranges, azm, major_axis, esquared);
}