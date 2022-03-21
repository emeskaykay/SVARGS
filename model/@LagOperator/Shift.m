% Shift
% Usage: Shift(shft)
% Shifts the lags in the Phi part of the operator given by:
%		Phi(L)(Y(t)) = A1.y(t-1) + A2.y(t-2) + ... + Ap.y(t-p)
% by a scalar amount s, ie:
%		Shift(Phi(L)(Y(t)), s) = Phi(L)(Y(t-s))
%
% A positive shift increases the lags, and vice versa.
% The maximum negative shift is Op.Lags(1)-1.
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
function Op = Shift(Op, shft)

	assert(shft > -Op.lags(1));
	Op.lags = Op.lags + shft;

	% Old Cache is no longer valid
	Op = CacheNew(Op);
	
end
%%
