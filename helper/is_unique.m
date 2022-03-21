% is_unique
%
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
function bu = is_unique(M, tol)
	if (nargin < 2) || isempty(tol), tol = 0; end
	assert(isnumeric(M));
	
	if numel(M) < 2
		bu = true;
		
	else
		M = vecc(M);
		dx = abs(diff(sort(M)));
		
		% By default, sort() places NaN's at the end for ascending sort.
		if isnan(dx(end))
			dx = dx(1:end-1);
		end
		
		bu = all(dx > tol);
	end
	
end
%%
