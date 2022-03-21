% rand_disk
%
% Generates random complex number(s) within the intersection of
% a specified annulus and a specified sector.
%
% rlims:		An (nx2) matrix of bounds on the absolute value (r)
%				of each complex number.
%
% ttalims:	An (nx2) matrix of bounds on the phase (tta). m must
%				be be less than or equal to n. If m is strictly less
%				than n, the remaining bounds on the phases are taken
%				to be [0, 2pi]. Any bound specified as [0, 0] will
%				generate a real number (whose sign is randomly chosen
%				to be positive or negative).
%
% Eigs:		Return value is a column vector of n complex numbers.

% The complex numbers are generated uniformly in (x,y) plane, over
% the area enclosed by the bounds.
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
function c = rand_disk(rlims, ttalims)
	if (nargin < 2), ttalims = []; end
	if (nargin < 1), rlims = [0, 1]; end
	
	n = size(rlims, 1);
	m = size(ttalims, 1);
	assert( numel(rlims) == 2*n );
	assert( numel(ttalims) == 2*m );
	ttalims = [ ttalims; repmat([0, pi], n-m, 1) ];

	% We need to take the square since we want to generate the
	% eigenvalues uniformly over the area in the (x,y) plane
	% (rather than in the (r,tta) plane). (ie: the number of
	% eigenvalues having length r should be proportional to r).
	rlimsqr = rlims.^2;
	rsqr = rlimsqr(:,1) + ( rlimsqr(:,2) - rlimsqr(:,1) ).*rand(n,1);
	r = sqrt(rsqr);
	
	tta = ttalims(:,1) + ( ttalims(:,2) - ttalims(:,1) ).*rand(n,1);
	
	signs = 2*(rand(n, 1) < 0.5) - 1;
	c = signs.*r.*complex(cos(tta), sin(tta));
	
end
%%
