% enforce_modifier
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
function enforce_modifier(n_out, n_min)

	if (nargin < 2) || isempty(n_min), n_min = 1; end
	
	if n_out < n_min
		try
			if n_min == 1
				ok = false;
				msg = 'Function is a modifier of a value class. Call it with an output argument. ';
				
			elseif n_min > 1
				ok = false;
				msg = ['Function requires at least ', num2str(n_min),  ' output arguments. '];
				
			else
				ok = true;
			end
			
			assert(ok, msg);
			
		catch ex
			throwAsCaller(ex);
		end
	end
	
end
%%
