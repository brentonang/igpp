% Function Cpath_MetProfiles		Reads profile data along a path
% all the data is in hdf5 format
%
% NOTE: this version does not do interpolation
%     [topoData,ifWater,olat,olon,orange] = Cpath_MetProfiles(path_lat,path_lon,ranges);
%   input:
%         path_lat - vector of latitudes
%         path_lon - vector of longitudes
%		  ranges - ranges from source
%   outpu
%         topoData - vector of topo data along the path
%         ifWater - vector gives 1 for land, 0 for water
%		  olat,olon - vectors where topo data is obtained (because no interpolation)
%		  orange - output ranges (nearest profiles)

function [topoData,ifWater,olat,olon,orange] = Cpath_MetProfiles(path_lat,path_lon,ranges);

DatabasesDir = 'C:/Users/User/Downloads/Google Drive/Miscellaneous/IGPP/2015072118';
filename = [DatabasesDir '/wrfout_d04_2015-07-21_12.GRM_F'];
%filename = [DatabasesDir '/wrfout_d03_2015-07-21_12.GRM_F'];
%filename = [DatabasesDir '/wrfout_d02_2015-07-21_12.GRM_F'];

% explore what is in the file
%dum = h5info(filename)
% see the results of this in diary_hdf5
%diary
%h5disp(filename)
%diary off

deg2rad = pi/180;
% convert to longitudes from -180 to 180 for this database
path_lon=Csetminmax(path_lon,-180,180);			

% read in grids of latitude, longitude, water (yes/no), terrain heights
% 
gridLat = h5read(filename,'/XLAT');				% latitude
gridLat_units = h5readatt(filename,'/XLAT','units');
gridLon = h5read(filename,'/XLONG');			% longitude
gridLon_units = h5readatt(filename,'/XLONG','units');

gridWater = h5read(filename,'/LANDMASK');   % LAND MASK (1 FOR LAND, 0 FOR WATER)
gridElev = h5read(filename,'/HGT');         % Terrain height
gridElev = h5read(filename,'/HGT');

% -------------------------------------------------------------------------------
% get topography data at nearest gridpoint for each point along the path
% only if the data are being loaded at a new grid point
% for each path point
%	a) find the distance from the point to each gridpoint
%	b) find the nearest gridpoint to the path point  
%	c) if it differs from the last gridpoint, copy the profile 

npath = length(path_lat); nprofs=0;
latpt = 100; lonpt = 400;
% -------------------------------------------------------------------------------
% for each point along the path
% -------------------------------------------------------------------------------
for jj=1:npath
% calculate the distance from the path point to all gridpoints in file
	coslat = cos(path_lat(jj)*deg2rad);
	xdist = coslat*(path_lon(jj)-gridLon);
	ydist = (path_lat(jj)-gridLat);
	dist = sqrt(xdist.*xdist+ydist.*ydist);
% find the gridpoint closest to the current path point
	dum = min(min(dist));
	[irow,jcol] = find(dist==dum);
	testlat = gridLat(irow,jcol); testlon = gridLon(irow,jcol);
	
% read the profile only if it is not at the same gridpoint as the last one
	if (testlat ~= latpt | testlon ~= lonpt)
		nprofs = nprofs+1;
		latpt = testlat; lonpt = testlon;	% reset point
		olat(nprofs) = testlat; olon(nprofs) = testlon;
		topoData(nprofs) = gridElev(irow,jcol);
        ifWater(nprofs) = gridWater(irow,jcol);
		
% now that I have a gridpoint, find the path point nearest to it
% (I need this extra step to get an approximate range value)		
		coslat = cos(testlat*deg2rad);
		xdist = coslat*(path_lon-testlon);
		ydist = (path_lat-testlat);
		dist = sqrt(xdist.*xdist+ydist.*ydist);	
		[dum,inr] = min(dist);	
		orange(nprofs) = ranges(inr);

%		[nprofs dum path_lat(jj) testlat path_lon(jj) testlon orange(nprofs)/1000]
	end
end 

%plot the entrie grid of latitudes and longitudes,
% with picked stations along the path
figure(11),clf,plot(gridLon,gridLat,'k.')
hold on,plot(olon,olat,'ro')
plot(path_lon(1),path_lat(1),'gp','MarkerSize',12,'MarkerFaceColor','g')
%keyboard