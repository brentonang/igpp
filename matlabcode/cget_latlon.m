%CGET_LATLON        Get latitudes and longitudes along a line
%
%       function [plat,plon]=cget_latlon(blat,blon,ranges,azm,major_axis);
%
% This is a translation to matlab of the USGS PROJ-4.4 geographic
% projection library to compute positions along a series of ranges
% along a given azimuth from a reference point.
% Uses ellipsoidal earth model set to Clark 1966 standard.
%
% Inputs:
%       blat,blon       lat,lon coordinates of start point
%                       in degrees (scalar)
%
%       ranges          range in nmi (scalar or vector)
%       azm             single-valued (same size as range) (degrees)
%       major_axes - skip for an elliptical earth, set=0 for spherical Earth
%
% Outputs:
%       plat,plon       lat,lons along path
%


function [plat,plon]=cget_latlon(blat,blon,ranges,azm,major_axis,esquared);

if nargin < 5
% Set the parameters for the "Clark 1966" ellipsoid model
   major_axis =  6378206.4;        % earth major axis in meters
   esquared   = 0.006768658;       % eccentricity squared
end

if major_axis == 0
   major_axis =  6378206.4;
   esquared   = 0;
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% now do the stuff that used to be in .c

len=length(ranges);
if length(azm)~=len
   disp('RANGES and AZM must be the same length')
   return
end
keyboard
deg2rad=pi/180;
rad2deg=180/pi;
halfpi=pi/2;
dtol=1e-12;
meritol=1e-10;

a=major_axis;S=ranges*1852;
phi1=blat*deg2rad;lam1=blon*deg2rad;al12=azm*deg2rad;
es=esquared;

if es ~= 0.0
   ellipse = 1;
else
   ellipse = 0;
end

if ellipse == 1
   onef= sqrt(1. - es);
   f  = 1. - onef;
   f4 = f/4.;
else
   onef = 1.;f=0;f4=0;
end

al12=cadjlon(al12);
is1=find(abs(al12) > halfpi);signS(is1)=1;
is0=find(abs(al12) <= halfpi);signS(is0)=0;
if ellipse == 1
   th1=atan(onef * tan(phi1));
else
   th1=phi1;
end
costh1 = cos(th1);sinth1 = sin(th1);
sina12 = sin(al12);
im1=find(abs(sina12) < meritol);merid(im1)=1;
im0=find(abs(sina12) >= meritol);merid(im0)=0;
for j=1:len
   if merid(j) == 1
      sina12(j) = 0.;
      if abs(al12(j)) < halfpi
         cosa12(j)=1;
      else
         cosa12(j)=-1;
      end
      M(j) = 0.;
   else 
      cosa12(j)= cos(al12(j));
      M(j)= costh1 * sina12(j);
   end
end
N = costh1 * cosa12;
if ellipse ==1
%   if merid(j) == 1
      c1(im1)= 0.; c2(im1)= f4; D(im1)= 1.-c2(im1); 
      D(im1) = D(im1).*D(im1); P(im1) = c2(im1)./D(im1);
%   else 
      c1(im0) = f*M(im0);c2(im0)= f4*(1. - M(im0).*M(im0));
      keyboard
      D(im0) = (1.-c2(im0)).*(1.-c2(im0)-c1(im0).*M(im0));
      P(im0) = (1.+.5*c1(im0).*M(im0)).*c2(im0)./D(im0);
%   end
end
for j=1:len
   if merid(j) == 1
      s1(j) = halfpi - th1;
   else 
      if abs(M(j)) >= 1
         s1(j) = 0;
      else
         s1(j) = acos(M(j));
      end
      s1(j)=  sinth1/sin(s1(j));
      if abs(s1(j)) >=1
         s1(j) = 0;
      else
         s1(j) = acos(s1(j));
      end
   end
end

if ellipse == 1
   d = S./(D * a); 
   d(is1) = -d(is1);
   u = 2.*(s1 - d);
   V = cos(u + d);
   sind = sin(d);
   X = c2.*c2.*sind.*cos(d).*(2. * V.*V - 1.);
   ds = d + X - 2. * P.*V.*(1. - 2. * P.*cos(u)).*sind;
   ss = s1 + s1 - ds;
else 
   ds = S/a;
   ds(is1) = - ds(is1);
end
cosds = cos(ds); sinds = sin(ds);
sinds(is1) = - sinds(is1);
al21 = N.*cosds - sinth1 * sinds;
for j=1:len
   if merid(j)==1
      phi2(j) = atan(tan(halfpi + s1(j)- ds(j))/onef);
      if al21(j) > 0.
         al21(j) = pi;
         if signS(j)==1; 
            de(j) = pi;
         else 
            phi2(j) = - phi2(j); de(j) = 0.;
         end              
      else 
         al21(j) = 0.;
         if signS(j)==1;
            phi2(j) = - phi2(j);de(j) = 0;
         else
            de(j) = pi;    
         end
      end          
   else 
      al21(j) = atan(M(j)/al21(j));
      if al21(j) > 0; al21(j) = al21(j) + pi; end
      if al12(j) < 0.; al21(j) = al21(j) - pi; end
      al21(j) = cadjlon(al21(j));
      if ellipse == 1
         denom = onef * M(j);
      else
         denom = M(j);
      end
      phi2(j)=atan(-(sinth1*cosds(j)+N(j)*sinds(j))*sin(al21(j))/denom);
      de(j)=atan2(sinds(j)*sina12(j),...
                 (costh1*cosds(j)-sinth1*sinds(j)*cosa12(j)));
      if ellipse == 1
         if signS(j)==1
            de(j)=de(j)+c1(j)*((1.-c2(j))*ds(j)+c2(j)*sinds(j)*cos(ss(j)));
         else
            de(j)=de(j)-c1(j)*((1.-c2(j))*ds(j)-c2(j)*sinds(j)*cos(ss(j)));
         end
      end
   end
end
keyboard
lam2 = cadjlon(lam1+de);
keyboard
lam1=lam1*rad2deg;
phi1 = phi1*rad2deg;
al12 = al12*rad2deg;
plon = lam2*rad2deg;
plat = phi2*rad2deg;

keyboard

