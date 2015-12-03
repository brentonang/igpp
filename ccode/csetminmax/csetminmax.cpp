#include <iostream>
#include <vector>
using namespace std;

vector<float> adjustToLimits (vector<float>& value, float lowerlim, float upperlim, float difference) {
	for (int i = 0; i < value.size(); i++) { // INSERT WHILE LOOP TO CHECK FOR VALUES GREATER THAN AZIMUTH
		if (value[i] >= upperlim) value[i] = value[i] - difference;
		else if (value[i] < lowerlim) value[i] = value[i] + difference;
	}
	return value;
}

void cSetminmax(vector<float>& value, float lowerlim, float upperlim) {
	float difference = upperlim - lowerlim;
	if (difference <= 0) cout << "Error: Upper limit must be higher than lower limit" << endl;
	adjustToLimits(value, lowerlim, upperlim, difference);
}

int main() {
	vector<float> myVector {1.0, 2.0, 3.0, 4.0, 5.0};
	cSetminmax(myVector, 2, 4);
	for (int i = 0; i < myVector.size(); i++) {
		cout << myVector[i] << endl;
	}
	return 0;
}