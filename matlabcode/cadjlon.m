%CADJLON        reduces argument to range from -pi to pi for single value, 
% use Csetminmax instead 
%
%       function [theta]=cadjlon(theta);
%

function [theta]=cadjlon(theta)

twopi=2*pi;

ind=find(theta>pi);
while isempty(ind)==0
	  theta(ind) = theta(ind)-twopi;
	  ind=find(theta>pi);
end

ind=find(theta< -pi);
while isempty(ind)==0
      theta(ind) = theta(ind)+twopi;
	  ind=find(theta< -pi);
end
