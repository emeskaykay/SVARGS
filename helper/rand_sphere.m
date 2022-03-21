% rand_sphere
% Generates a random variable that is uniformly ditributed
% inside a k-hyperspere between radii ri and ro.
%
% k:	Dimension of the hypersphere.
% ro: Outer radius.
% ri: Inner radius.
%
% * If ri == ro, then a random point on the surface of radius ri is returned.
% * Default values: ro = 1, ri = 0.
% * ro and ri can be vectors (of the same length n), in which case n points
%		are returned, the j-th point lying between ri(j) and ro(j).
%
% Notes:
%	The implementation uses k independent normal standard random variables.
%	Since the joint distribution of these k variables is spherically symmetric,
%	one only has to normalize by the appropriate radius to get a spherically
%	symmetric distribution.
%
%	Other implementations are also possible. For example, one could generate
%	k uniform variables (x1, x2, ..., xk) on [0, 1] and then reject any points
%	for which x1^2 + x2^2 + ... xk^2 > 1. However this method will be
%	inefficient for large values of k.
%
% For other implementations, see for example:
%		"Choosing a point from the surface of a sphere, by George Marsaglia.
%
% For the special case of a 2-sphere, the one can choose a uniform random angle
% between 0 and 2pi. However this method does not carry over easily for the case
% of k >= 3.
%
%%
% Copyright 2021, Sandeep K Mody, Govindan Rangarajan
% This file is part of the SVARGS program.
% SVARGS is free software: you can redistribute it and/or modify it under the terms
% of the GNU Lesser General Public License as published by the Free Software Foundation,
% either version 2 of the License, or (at your option) any later version.
%
% SVARGS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
% PURPOSE. See the GNU Lesser General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public License along with
% SVARGS. If not, see <https://www.gnu.org/licenses/>.
% 
%%
function v = rand_sphere(k, ro, ri)
	if (nargin < 3), ri = 0; end
	if (nargin < 2), ro = 1; end
	
	n = length(ri);
	assert(numel(ri) == length(ri));
	assert(numel(ro) == n);
	
	if (size(ri, 1) > size(ri, 2))
		ri = ri';
	end
	
	if (size(ro, 1) > size(ro, 2))
		ro = ro';
	end
	
	si = ri.^k;
	so = ro.^k;
	s = si + (so - si).*rand(1,n);
	r = s.^(1/k);
	
	x = randn(k, n);
	l = sqrt(sum(x.*x, 1));
	v = zeros(k, n);
	for i = 1 : k
		x(i, :) = x(i, :)./l;
		v(i,:) = r.*x(i,:);
	end
	
end

%%
% TEST
% v = rand_sphere(3);
% plot3(v(1), v(2), v(3), '*');
% hold on;
% for i = 1 : 1000
% 	v = rand_sphere(3);
% 	plot3(v(1), v(2), v(3), '*');
% end
%
%%
