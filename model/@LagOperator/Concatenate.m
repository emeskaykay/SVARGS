% Concatenate
% Concatenate LagOperator objects
%
% dim:
%	The dimension along which to concatenate.
%		1 -> diagonal concatenation along variables, 2 -> along lags.
%
% ops:
%	Vector of lag operator objects.
%
% Returns:
%	Concatenated LagOperator object.
%
%% TODO: Implement for dim > 1
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
function Obj = Concatenate(dim, ops)

	assert(nargin > 1);
	assert( isnumeric(dim) );
	assert(dim == 1, 'Concatenate not yet implemented for dim > 1.');
	assert( isa(ops, 'LagOperator' ) );
	
	extlags = ops(1).Lags();
	nmods = numel(ops);
	for i = 2 : nmods
		extlags = union(extlags, ops(i).Lags());
	end

	m = length(extlags);
	kvec = arrayfun(@NumVars, ops);
	ksum = sum(kvec);
	Cdiag = zeros(ksum, ksum, m);
	Jdiag = zeros(ksum, ksum);
	kcumsum = 0;
	for i = 1 : nmods
		r = kcumsum+1 : kcumsum + kvec(i);
		Jdiag(r,r) = ops(i).Jmat;
		for j = 1 : m
			Cdiag(r,r,j) = Coeffs(ops(i), extlags(j));
		end
		kcumsum = kcumsum + kvec(i);
	end
	Obj = LagOperator(extlags, Cdiag, 'J', Jdiag);
	
end
%%
