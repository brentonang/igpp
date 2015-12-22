srclat = 41.131; srclon = 360-112.8965;		% source location
rinc = 100;									% range increment (m)
azimuth = 70;								% azimuth cw from N (degrees)
range = 160000;								% maximum range (m)
maxrange = 130000;

%path(path,'C:/Users/User/Downloads/Google Drive/Miscellaneous/IGPP')
%path(path,'/Users/cdh/Catherine/CH2MHill/topography');

% .......................................................................
% 1) get vector of ranges along the path (maxrange <=130 km)
% .......................................................................
max_range = min(maxrange,range);
ranges = [0:rinc:max_range+rinc];

% .......................................................................
% 2) get latitudes and longitudes along the path (must be eastward) .....
% .......................................................................
% correct azimuth to an eastward value
% see footnote 1	to convert to C++
if azimuth < 0 | azimuth >180
	disp('ERROR: azimuth must be Eastward')
	return
end
azm = azimuth*ones(size(ranges));
[plat,plon]=cget_latlon(srclat,srclon,ranges/1852,azm); % convert to nmi

plon=Csetminmax(plon,0,360);			% convert to positive longitudes

keyboard
% do a quick comparison of path coords on geodesic surface vs. flattish Earth

flatlon = srclon+sind(azimuth)*ranges/(cosd(srclat)*111111);
flatlat = srclat+cosd(azimuth)*ranges/111111;
figure(11),clf,
plot(plon,plat,'k'),hold on, plot(flatlon,flatlat,'r')

keyboard
% ......................................................................
% 3) get the topography along the path
% ......................................................................

[ptopo] = Cpath_topoSRTM1(plat,plon);

y=[ranges; ptopo];
fname= ['topo',int2str(azimuth),'.in']
fid = fopen(fname,'w');
fprintf(fid,'%11.2f  %12.1f\n',y);
fclose(fid);

figure(1),clf,plot(y(1,:)/1000,y(2,:))
title(['Topography along ',int2str(azimuth),' degrees azimuth'])
xlim([0 y(1,end)/1000]),grid
xlabel('range (km)')
ylabel('Altitude (m)')


return
		
% ptopo gives depths to bottom of the lake
% it needs to be amended to account for water in lake
% as is, it assumes that the lake is dry
% ...........................................................

% compare to values from sparser global database
%[topo_data,olat,olon] = Cpath_topo30(plat,plon);

% *********************************************************************

[CH2Mtopo,ifWater,olat,olon,orange] = Cpath_topoCH2M(plat,plon,ranges);


figure(1),clf,axes('position',[0.1 0.3 0.8 0.6])
plot(ranges/1000,topo_data,'k'),grid
hold on, plot(ranges/1000,ptopo,'r')
hold on, plot(orange/1000,CH2Mtopo,'b')
ind = find(ifWater==0);
plot(orange(ind)/1000,CH2Mtopo(ind),'bo')
legend('topo30','SRTM1','CH2MHill','water','Location','NorthWest')
title(['Topography along azimuth ',num2str(azimuth),' degrees'])
ylabel('Altitude (m)')
xlim([0 maxrange/1000])

axes('position',[0.1 0.1 0.8 0.15])
plot(orange/1000,ifWater,'-bo'),grid
ylabel('water/land')
xlabel('Range (km)')
xlim([0 maxrange/1000])

print -djpeg99  compareTopos.jpg