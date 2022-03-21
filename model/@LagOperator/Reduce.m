% Op = Reduce(Op)
%
%
%% TODO: check if any cached variables are affected.
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
function Op = Reduce(Op)
	enforce_modifier(nargout);
	
	if Op.has_lag_zero
		[Op.Phi, ~, Op.Jmat] = Op.Coeffs(Op.lags, 'reduced', true);
	end
end
%%
