#include "csetminmax.hpp"

#include <iostream>

using namespace std;

vector<double> adjustToLimits(vector<double> value, double lowerlim, double upperlim, double difference) {
	for(int i = 0; i < value.size(); i++) {
		while(value[i] >= upperlim || value[i] < lowerlim) {
			value[i] = (value[i] >= upperlim) ? value[i] - 2*difference : value[i] + 2*difference;
		}
	}
	return value;
}

double adjustToLimits(double value, double lowerlim, double upperlim, double difference) {
	while(value >= upperlim || value < lowerlim) {
		value = (value >= upperlim) ? value - difference : value + difference;
	}
	return value;
}

vector<double> cSetMinMax(vector<double> value, double lowerlim, double upperlim) {
	double difference = upperlim - lowerlim;
	if (difference <= 0) {
		cout << "Error: Upper limit must be higher than lower limit" << endl;
		vector<double> zeros(value.size());
		return zeros;
	}
	return adjustToLimits(value, lowerlim, upperlim, difference);
}

double cSetMinMax(double value, double lowerlim, double upperlim) {
	double difference = upperlim - lowerlim;
	if (difference <= 0) {
		cout << "Error: Upper limit must be higher than lower limit" << endl;
		return 0;
	}
	return adjustToLimits(value, lowerlim, upperlim, difference);
}