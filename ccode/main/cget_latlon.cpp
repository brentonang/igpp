#include "csetminmax.hpp"
#include "cget_latlon.hpp"

#include <algorithm>
#include <iostream>
#include <utility>
#include <fstream>

using namespace std;

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

pair<vector<double>, vector<double>> cGetLatLon(double blat, double blon, vector<double> ranges, vector<double> azm, double major_axis = 6378206.4, double esquared = 0.006768658) {
	double costh1, denom, ellipse, f, f4, lam1, onef, phi1, th1, sinth1;
	int i, j, len = ranges.size();
	pair<vector<double>, vector<double>> returnValues;
	vector<double> al12(1302), al21(1302), c1(1302), c2(1302), cosa12(1302), cosds(1302), D(1302), d(1302), de(1302), ds(1302), lam2(1302), M(1302), merid(1302), N(1302), P(1302), phi2(1302), plon(1302), plat(1302), S(1302), s1(1302), ss(1302), signS(1302), sina12(1302), sind(1302), sinds(1302), u(1302), V(1302), X(1302);

	if (major_axis == 0) {
		major_axis = 6378206.4;
		esquared = 0;
	} 

	if (azm.size() != len) {
		cout << "Ranges and AZM must be the same length" << endl;
		vector<double> zeros1(len);
		vector<double> zeros2(len);
		returnValues.first = zeros1;
		returnValues.second = zeros2;
		return returnValues;
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

	al12 = cSetMinMax(al12, -180, 180);
	for (i = 0; i < al12.size(); i ++) { 
		sina12[i] = sin(al12[i]);
		cosa12[i] = cos(al12[i]);			
		signS[i] = (abs(al12[i]) > halfpi) ? 1 : 0;
	}

	th1 = phi1;
	if (ellipse == 1) th1 = atan(onef * tan(phi1));
	costh1 = cos(th1);
	sinth1 = sin(th1);

	for (i = 0; i < len; i++) {
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
				c2[i] = f4 * (1 - M[i] * M[i]);
				D[i] = (1 - c2[i]) * (1 - c2[i] - c1[i] * M[i]);
				P[i] = (1 + 5 * c1[i] * M[i]) * (c2[i]/D[i]);
			}
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
	for (j = 0; j < len; j++) {
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
			if (al12[j] < 0) al21[j] -= pi; //al12 or al21?
			al21[j] = cSetMinMax(al21[j], -pi, pi);
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
	for (i = 0; i < len; i++) lam2[i] = cSetMinMax(de[i], -pi, pi);
	lam1 *= rad2deg;
	phi1 *= rad2deg;
	al12 = elementWiseMultiplication(al12, rad2deg);
	plon = elementWiseMultiplication(lam2, rad2deg);
	plat = elementWiseMultiplication(phi2, rad2deg);

	returnValues.first = plat;
	returnValues.second = plon;
	return returnValues;
}

int main() {
	int i = 0, j = 0;
	ofstream myFile;
	double max_range, srclat = 41.131, srclon = 360-112.8965, rinc = 100, azimuth = 70, range = 160000, maxrange = 130000;
	vector<double> tempValues(1302), ranges(1302), plat(1302), plon(1302); 
	pair<vector<double>, vector<double>> myPair;
	max_range = min(maxrange, range);
   while(j <= max_range + rinc) {
      ranges[i] += j;
      j += rinc;
      i++;
   } 

   vector<double> azm(ranges.size(), 1 * azimuth);

   tempValues = ranges;
	for (i = 0; i < ranges.size(); i++) {
		tempValues[i] = ranges[i] / 1852;
	}
	myPair = cGetLatLon(srclat, srclon, tempValues, azm, 6378206.4, 0.006768658);
	plat = myPair.first;
	plon = myPair.second;
	myFile.open("cgetlatlonoutput.txt");
	myFile << "PLAT VALUES ARE: " << endl;
	for(vector<double>::iterator it = plat.begin(); it != plat.end(); it++) myFile << *it << " ";
	myFile << endl << endl << "PLON VALUES ARE: " << endl;
	for(vector<double>::iterator it = plon.begin(); it != plon.end(); it++) myFile << *it << " ";
	myFile.close();
	return 0;
}