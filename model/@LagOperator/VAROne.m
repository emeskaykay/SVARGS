% VAROne
% Returns the corresponding AR(1) Operator from the AR(p) operator.
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
function Op1 = VAROne(Op)
	assert( isa(Op, 'LagOperator') );
	
	p = Op.Order();
	Phi1 = BlockCompanion(Op.Coeffs(1:p, 'reduced', true));
	if isempty(Phi1)
		lags = [];
	else
		lags = 1;
	end
	Op1 = LagOperator(lags, Phi1);
end
%%
