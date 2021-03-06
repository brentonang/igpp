#include "csetminmax.hpp"
#include "cget_latlon.hpp"
#include "Cpath_topoSRTM1.hpp"
#include "CompareTopo_AlongPath.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <tuple>
#include <vector>
#include <utility>

using namespace std;

int main(int argc, char* argv[]) {
   ofstream myFile;
   double azimuth, maxrange, max_range, range = 160000, rinc, scalar1, scalar2, srclat = 41.131, srclon = 360 - 112.8965;
   int i = 0, j = 0;

   if(argc != 4) {
      cout << "Invalid number of command line arguments" << endl;
      return -1;
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
   
   int vector_size = floor((maxrange + rinc) / rinc) + 1;
   vector<double> flatlat(vector_size), flatlon(vector_size), olat(vector_size), olon(vector_size), plat(vector_size), plon(vector_size), ranges(vector_size), tempRanges(vector_size), vector1(vector_size), vector2(vector_size);
   pair<vector<double>, vector<double>> tempValues;
   pair<vector<double>, vector<double>> returnPair;
   // tuple<vector<double>, vector<double>, vector<double>> returnTuple;

   max_range = min(maxrange, range);
   while(j <= max_range + rinc) {
      ranges[i] += j;
      j += rinc;
      i++;
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
   for (i = 0; i < plon.size(); i++) plon[i] = cSetMinMax(plon[i], 0, 360);
   for (i = 0; i < ranges.size(); i++) flatlon[i] = srclon + sin(azimuth * deg2rad) * ranges[i] / (cos(srclat * deg2rad) * 111111);
   for (i = 0; i < ranges.size(); i++) flatlat[i] = srclat + cos(azimuth * deg2rad) * ranges[i] / 111111;

   returnPair = Cpath_topoSRTM1(plat, plon);
   olat = returnPair.first;
   olon = returnPair.second;
   myFile.open("igppoutput.txt");
   myFile << "OLAT [" << olat.size() << "] VALUES ARE: " << "\n";
   for(vector<double>::iterator it = olat.begin(); it != olat.end(); it++) myFile << std::setprecision(5) << *it << "\t\t";
   myFile << "\n" << "\n" << "OLON [" << olon.size() << "] VALUES ARE: " << "\n";
   for(vector<double>::iterator it = olon.begin(); it != olon.end(); it++) myFile << std::setprecision(5) << *it << "\t\t";
   myFile << "\n" << "\n" << "PLAT [" << plat.size() << "] VALUES ARE: " << "\n";
   for(vector<double>::iterator it = plat.begin(); it != plat.end(); it++) myFile << std::setprecision(5) << *it << "\t\t";
   myFile << "\n" << "\n" << "PLON [" << plon.size() << "] VALUES ARE: " << "\n";
   for(vector<double>::iterator it = plon.begin(); it != plon.end(); it++) myFile << std::setprecision(5) << *it << "\t\t";
   myFile.close();
   cout << "Calculation finished. Output values are written to igppoutput.txt" << endl;
   return 0;
}
