% Function CSETMINMAX	resets value to within given limits
%	useful for ensuring that longitude is in correct range
%      [value] = Csetminmax(value,lowlim,uplim)
%	
%	input:
%		value - value to be limited
%		lolim = lower limit
%		uplimit = upper limit
%	output:
%		value - new value (between limits)
%
%	EXAMPLE:
%	a) force longitude 0<=longitude<360
%		longitude = -180;
%		[longitude] = Csetminmax(longitude,0,360)
%	b) force longitude -pi<=longitude<pi
%		longitude = 2*pi;
%		[longitude] = Csetminmax(longitude,-pi,pi)

function	[value] = Csetminmax(value,lowlim,uplim);

diff=uplim-lowlim;
if diff <=0
	disp('ERROR: upper limit must be higher than lower limit')
	return
end

ind = find(value>=uplim);
while isempty(ind)==0
	value(ind) = value(ind)-diff;
    ind = find(value>=uplim);
end

ind = find(value<lowlim);
while isempty(ind)==0
	value(ind) = value(ind)+diff;
	ind = find(value<lowlim);
end

return