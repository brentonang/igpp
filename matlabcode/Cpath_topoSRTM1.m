% Function Cpath_topoSRTM1		Reads topography data along a path
% with given coordinates from the N39W113_4X3 file
% eventually may want to allow different topo files
%
% NOTE: this version does not do interpolation
%     [path_data,olat,olon] = Cpath_topoSRTM1(path_lat,path_lon);
%   input:
%         path_lat - vector of latitudes
%         path_lon - vector of longitudes
%   output:
%         path_data - topography from SRTM1
%		  olat,olon - vectors where topo data is obtained (because no interpolation)
% 
function  [path_data,olat,olon] = Cpath_topoSRTM1(path_lat,path_lon);

rad2deg=180/pi; deg2rad=pi/180;
% databases directory location changed from global to local variable
DatabasesDir = '/Users/cdh/Catherine/CH2MHill/topography';
fname = [DatabasesDir '/N39W113_4X3.hgt'];
% 
% NOTE: the data in this file are ordered as
% first 14401 values are from S to N along the first longitude (247E)
% etc
% last 14401 values are from S to N along the last longitude (250E)

% Set up parameters defining the N39W113_4X3.hgt file 
% parameters are hard-wired to this particular file
db_res         = 1/3600;		% 1 second resolution
hdb_res=db_res/2;
db_loc         = [39 43 247 250];	% region covered by topography file
db_size        = [14401 10801];		% # of points in latitude & longitude
nbytes_per_lon = db_size(1)*2;	% 2-byte integers

% Read in the data directly
	
	% Make sure the longitude data goes from 0 to 360
	path_lon = Csetminmax(path_lon,0,360);

	% Calculate the lat/lon indices into the matrix
 	ilat = round((path_lat - db_loc(1))/db_res);
	ilon = round((path_lon  - db_loc(3))/db_res); 

	% Calculate the offset to each point
	offset = ilon*nbytes_per_lon + ilat*2;

	% Open the data file  (drop ieee-be)
	fid = fopen(fname, 'r');
	if (fid < 0)
	  errordlg(['Could not open topo file: ' DatabasesDir '/N39W113_4X3.hgt'],'Error');
	end

	% Go ahead and read the database
	path_data = zeros(size(path_lon));
	for i = 1:length(path_lon);
		status = fseek(fid, offset(i), 'bof');
		path_data(i)=fread(fid,[1],'integer*2');
	end

fclose(fid);							% close the file

glon = [1:db_size(2)]*db_res-hdb_res;
glat = 90-[1:db_size(1)]*db_res+hdb_res;
olon=glon(min(db_size(2),ilon+1))-hdb_res;
olat=glat(min(db_size(1),ilat+1))+hdb_res;

