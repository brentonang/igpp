#include "Cpath_topoSRTM1.hpp"
#include "csetminmax.hpp"
#include "cget_latlon.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>  

using namespace std;

vector<int> returnMin(vector<int> firstVector, int secondValue) {
   vector<int> returnVector(firstVector.size());
   for (int i = 0; i < firstVector.size(); i++) returnVector[i] = (firstVector[i] < secondValue) ? firstVector[i] : secondValue;
   return returnVector;
}

// tuple<vector<double>, vector<double>, vector<double>> Cpath_topoSRTM1(vector<double> path_lat, vector<double> path_lon) {
pair<vector<double>, vector<double>> Cpath_topoSRTM1(vector<double> path_lat, vector<double> path_lon) {
   int i = 0, j = 0, status, len = path_lat.size();
   long fileSize;
   double nbytes_per_lon;
   vector<int> db_loc{39, 43, 247, 250}, db_size{14401, 10801};
   vector<double> glon(db_size[1]), glat(db_size[0]), olon(len), olat(len), path_data(path_lon.size());
   vector<int> ilat(len), ilon(len), ilatplus1(len), ilonplus1(len), minLat(len), minLon(len), offset(len);
   tuple<vector<double>, vector<double>, vector<double>> returnTuple;
   tuple<vector<double>, vector<double>, vector<double>> zeroTuple;
   pair<vector<double>, vector<double>> returnPair;

   nbytes_per_lon = 2.0 * db_size[0];

   path_lon = cSetMinMax(path_lon, 0, 360);
   for (i = 0; i < path_lon.size(); i++) {
     ilat[i] = round((path_lat[i] - db_loc[0])/db_res);
     ilon[i] = round((path_lon[i] - db_loc[2])/db_res);
     offset[i] = ilon[i] * nbytes_per_lon + ilat[i] * 2;
   }

   ifstream file("N39W113_4X3.hgt", ios::in|ios::binary);
   if(!file) {
     cout << "Error opening file!" << endl;
     vector<double> zeros1(path_lat.size());
     vector<double> zeros2(path_lat.size());
     // vector<double> zeros3(path_lat.size());
     // std::get<0>(zeroTuple) = zeros1;
     // std::get<1>(zeroTuple) = zeros2;
     // std::get<2>(zeroTuple) = zeros3;
     // return zeroTuple;
     returnPair.first = zeros1;
     returnPair.second = zeros2;
     return returnPair;
   }

   const int SRTM_VERSION = 1201;
   int height[SRTM_VERSION][SRTM_VERSION];
   unsigned char memblock[2];
   for (i = 0; i < SRTM_VERSION; i++) {
    for (j = 0; j < SRTM_VERSION; j++) {
      height[i][j] = (memblock[0] << 8) | memblock[1];
    }
   }

   const int row = 500;
   const int column = 1000;
   size_t tempoffset = sizeof(memblock) * ((row * SRTM_VERSION) + column);
   file.seekg (tempoffset, ios::beg);
   file.read (reinterpret_cast<char*>(memblock), sizeof(memblock));
   short single_value = (memblock[0] << 8) | memblock[1];
   cout << "values at " << row << ", " << col << ": " << endl;
   cout << "  height array: " << height[row][col] << ", file: " << single_value << endl;

   // for(int i = 0; i < ilat.size(); i++) {
   //    for(int j = 0; j < ilon.size(); j++) {
   //       size_t bufferoffset = sizeof(buffer) * ((i * ilat.size()) + j);
   //       file.seekg(bufferoffset, ios::beg);
   //       file.read( reinterpret_cast<char*>(buffer), sizeof(buffer) );
   //       path_data[i * ilat.size() + j] = height[ilat[i]][ilon[j]];
   //    }
   // }

   for (i = 0; i < db_size[1]; i++) glon[i] = (i + 1) * db_res - hdb_res;
   for (int j = 0; j < db_size[0]; j++) glat[j] = 90 - (j + 1) * db_res + hdb_res;
   for (int k = 0; k < ilat.size(); k++) ilatplus1[k] = ilat[k] + 1;
   for (int l = 0; l < ilon.size(); l++) ilonplus1[l] = ilon[l] + 1;
   minLon = returnMin(ilonplus1, db_size[1]);
   minLat = returnMin(ilatplus1, db_size[0]);
   for (int m = 0; m < minLon.size(); m++) olon[m] = glon[minLon[m]] - hdb_res;
   for (int n = 0; n < minLat.size(); n++) olat[n] = glat[minLat[n]] + hdb_res;
   // std::get<0>(returnTuple) = path_data;
   // std::get<1>(returnTuple) = olat;
   // std::get<2>(returnTuple) = olon;
   // return returnTuple;
   returnPair.first = olat;
   returnPair.second = olon;
   return returnPair;
}

int main() {
   int i = 0, j = 0;
   ofstream myFile;
   double azimuth = 70, maxrange = 130000, max_range, range = 160000, rinc = 100, srclat = 41.131, srclon = 360 - 112.8965;
   vector<double> flatlat(1302), flatlon(1302), olat(1302), olon(1302), path_data(1302), plat(1302), plon(1302), ranges(1302), tempRanges(1302);
   pair<vector<double>, vector<double>> tempValues;
   tuple<vector<double>, vector<double>, vector<double>> returnTuple;
   pair<vector<double>, vector<double>> returnPair;

   max_range = min(maxrange, range);
   while(j <= max_range + rinc) {
     ranges[i] += j;
     j += rinc;
     i++;
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
   // returnTuple = Cpath_topoSRTM1(plat, plon);
   // path_data = std::get<0>(returnTuple);
   // olat = std::get<1>(returnTuple);
   // olon = std::get<2>(returnTuple);
   cout << "Received path_data, olat, and olon values" << endl;
   cout << "Creating file" << endl;
   myFile.open("cpathoutput.txt");
   cout << "Writing olat values, please wait" << endl;
   myFile << "OLAT VALUES ARE: " << endl;
   for(vector<double>::iterator it = olat.begin(); it != olat.end(); it++) myFile << *it << " ";
   cout << "Writing olon values, please wait" << endl;
   myFile << endl << endl << "OLON VALUES ARE: " << endl;
   for(vector<double>::iterator it = olon.begin(); it != olon.end(); it++) myFile << *it << " ";
   // cout << "Writing path_data values, please wait" << endl;
   // myFile << endl << endl << "PATH_DATA VALUES ARE: " << endl;
   // for (vector<double>::iterator it = path_data.begin(); it != path_data.end(); i++) myFile << *it << " ";
   // myFile.close();
   return 0;
}
