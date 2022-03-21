% CheckBinaryOp__
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
function [left, right] = CheckBinaryOp__(left, right, opsym)

	ERR_SIZEMATCH = ['Operator (', opsym, '): RHS and LHS have incompatible sizes. '];
	ERR_TYPE = 'LHS of (*) must be LagOperator or a numeric matrix or scalar. ';

	assert(isa(right, 'LagOperator'), ['Binary (', opsym, ') RHS must be LagOperator. ']);
	
	if (opsym == '*')
		if isa(left, 'LagOperator')
			assert( left.NumVars() == right.NumVars(), ERR_SIZEMATCH );
		elseif isnumeric(left)
			if ismatrix(left) && ~isscalar(left)
				assert( all(size(left) == right.NumVars()), ERR_SIZEMATCH );
			elseif ~isscalar(left)
				assert(false, ERR_TYPE);
			end
		else
			assert(false, ERR_TYPE);
		end
		
	else
		assert(isa(left, 'LagOperator'), ['Binary (', opsym, ') LHS must be LagOperator. ']);
		assert( left.NumVars() == right.NumVars(), ERR_SIZEMATCH );
	end
	
end
%%
