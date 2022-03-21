% is_integral
% Checks if all elements of the supplied array are finite integers, optionally
% within the specified limits and optionally of a given integer type.
%
% A:
%	Numeric array.
%
% lims [optional]:
%	A pair [a,b] specifying the minimum and maximum permissible values in A.
%
% classname [optional]:
%	String specifying the type of A. For example: 'double', 'int32' etc.
%
% OUTPUT:
% yes:
%	True if all elements of A are integers and satisfy the specified conditons.
%
% B:
%	Elementwise true or false.
%
% NOTE: A must be numeric or an error is thrown.
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
function [yes, B] = is_integral(A, lims, classname)
	assert(isnumeric(A));
	
	if (nargin < 3), classname = []; end
	if (nargin < 2) || isempty(lims), lims = [-Inf, Inf]; end
	
	if ~isempty(classname) && ~isa(A, classname)
		B = false(size(A));
		yes = false;
		return;
	end
	
	B = isfinite(A);
	
	if ~isinf(lims(1))
		B = B & (A >= lims(1));
	end
	
	if ~isinf(lims(2))
		B = B & (A <= lims(2));
	end
	
	B = B & all(rem(A,1) == 0);
	yes = all(vecc(B));
end
%%
