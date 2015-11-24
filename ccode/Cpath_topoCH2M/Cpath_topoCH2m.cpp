#include <iostream>
#include <ifstream>
#include <math.h>
#include <vector>

// #include "Cpath_topoCH2m.h"
using namespace std;

const double PI = 3.141592653589793;
// const double PI = 4*atan(1); // equates only to 3.14159

/****************************************************************************************************
* Function Cpath_MetProfiles	Reads profile data along a path
* all the data is in hdf5 format
*
* NOTE: this version does not do interpolation
*	[topoData, ifWater, olat, olon, orange] = Cpath_MetProfiles(path_lat,path_lon,ranges);
*	input:
*		path_lat - vector of latitudes
*		path_lon - vector of longitudes
*		ranges - ranges from source
*	output:
*		topoData - vector of topo data along the path
*		ifWater - vector gives 1 for land, 0 for water
*		olat, olon - vectors where topo data is obtained (because no interpolation)
*		orange - output ranges (nearest profiles)
**************************************************************************************************/

int Cpath_MetProfiles(vector<double> path_lat, vector<double> path_lon, vector<double> ranges) {
	ifstream fileName;
	fileName.open("wrfout_d04)2015-07-21_12");

	double deg2rad = PI/180;
	// convert to longitudes from -180 to 180 for this database
	path_lon = Csetminmax(path_lon,-180,180);

	/* read in grids of latitude, longitude, water(yes/no), terrain heights
	 * Need to research how to read these files in C++
	 *gridLat = h5read(filename,'/XLAT');				% latitude
	 *gridLat_units = h5readatt(filename,'/XLAT','units');
	 *gridLon = h5read(filename,'/XLONG');			% longitude
	 *gridLon_units = h5readatt(filename,'/XLONG','units');

	 *gridWater = h5read(filename,'/LANDMASK');   % LAND MASK (1 FOR LAND, 0 FOR WATER)
	 *gridElev = h5read(filename,'/HGT');         % Terrain height
	 *gridElev = h5read(filename,'/HGT');
	 */

	 /*********************************************************************************************
	  *get topography data at nearest gridpoint for each point along the path
	  *only if the data are being loaded at a new grid point
	  *for each path point
	  *a) find the distance from the point to each gridpoint
	  *b) find the nearest gridpoint to the path point  
	  *c) if it differs from the last gridpoint, copy the profile 
	  ********************************************************************************************/

	int npath = length(path_lat);
	int nprofs = 0;
	int latpt = 100;
	int lonpt = 400;

	for (int i = 1; i < npath; i++) {
		double coslat = cos(path_lat(i)*deg2rad);
		double xdist = coslat*(path_lon(i)-gridLon);
		double ydist = (path_lat(i)-gridLat);
		// double dist = sqrt(xdist*xdist+ydist*ydist);
		// need to figure out element wise multiplication
		// current plan is std::transform(v1.begin(), v1.end(), v2.begin(), v.begin(), std::multiplies<double>());

		// find the gridpoint closes to the current path point
		int dum = min(min(dist));
		// NEED TO FIX THIS
		// [irow, jcol] = find(dist==dum)
		// double testlat = gridLat(irow, jcol);
		// double testlon = gridLon(irow, jcol);

		// read the profile only if it is not at the same gridpoint as the last one
		if (testLat != latpt || testLon != lonpt) {
			nprofs++;
			latpt = testlat;
			lonpt = testlon;
			vector<double> olat, olon;
			olat(nprofs) = testlat;
			olon(nprofs) = testlon;
			topoData(nprofs) = gridEleve(irow, jcol);
			ifWater(nrpofs) = gridWater(irow, jcol);

			costlat = cos(testlat*deg2rad);
			xdist = coslat*(path_lon-testlon);
			ydist = (path_lat-testlat);
			// dist = sqrt(xdist.*xdist+ydist.*ydist);
			// [dum,inr] = min(dist);
			// orange(nprofs) = ranges(inr)
		}
	}

	// plot the entire grid of latitudes and longitudes
	// with picked stations along the path
	// figure(11),clf,plot(gridLon,gridLat,'k.')
	// hold on,plot(olon,olat,'ro')
	// plot(path_lon(1),path_lat(1),'gp','MarkerSize',12,'MarkerFaceColor','g')
}
