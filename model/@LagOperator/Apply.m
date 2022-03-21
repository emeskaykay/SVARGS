% Z = Apply(Op, Y, itimerange)
% Uses the operator to propagate values of Y to each time in the specified range,
%
% Y:
%	A (k x len) data vector, where k is the number of variables and 1:len is the
%	set of time points.
%
% iknownvars [Optional, default: empty]:
%	Variables indices of the variables that are already known at the result time.
%	A non empty vector for this argument can be used for performing "Nowcasts". The
%	indices in this vector must be a proper subset of (1:k) where k is the total
%	number of variables (number of rows in Y).
%
% itime_range [ Optional, default: [1,len] ]:
%	This specifies the range of times to which the operator should be propagated.
%	For each index t in the specified range, values of Y _prior_ to t as well as
%	the _known_ values of Y _at_ t (based on the indices supplied in iknownvars)
%	are used to propagate the (remaining) values of Y at t. Thus, t is the time at
%	which the values of Y are required, (t-1) is treated as lag 1, (t-2) as lag 2
%	and so on, upto the last lag p corresponding to index t-p. If t is less than or
%	equal to the order p, the data Y is padded at the start with p-t+1 zero columns.
%
%	The values in Y(ivars_unknow, end) (ivars_unknow == setdiff(1:k, iknownvars))
%	are ignored. They should be set to 0 by the caller.
%
% Z:
%	A (k x n) array, where n is the length of itime_range.
%	For a VAR model with this operator, the results are nowcasts or forecasts.
%	When iknownvars has at least one index, Z(:,i) can be regarded as the
%	nowcast at time i given Y(:,j) for j < i and given Y(iknownvars, i).
%	Otherwise Z(:, i) can be regarded as the 1-step forecast at time i given
%	Y(:,j) for j < i. Note that for nowcasts, the known variables in Y will
%	simply be copied into Z. Only the unkown variables will have values written.
% 
% See also: VARModel:Forecast().
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
function Z = Apply(Op, Y, iknownvars, itime_range)
	[k, n] = size(Y);
	assert( k == Op.NumVars(), 'LagOperator::Apply: Wrong number of variables. ' );
	
	if (nargin < 4) || isempty(itime_range), itime_range = [1, n]; end
	assert(isnumeric(itime_range) && (numel(itime_range) <= 2));
	assert( all(itime_range >= 0) && all(itime_range <= n) );
	if isscalar(itime_range), itime_range = [itime_range, itime_range]; end
	assert(itime_range(1) <= itime_range(2));
	
	if (nargin < 3) || isempty(iknownvars), iknownvars = zeros(1,0); end
	assert(isnumeric(iknownvars) && isvector(iknownvars));
	s = numel(iknownvars);
	assert(s < k);
	iknownvars = vecc(iknownvars);
	
	first = itime_range(1);
	last = itime_range(2);
	Z = Y(:, first:last);
	
	if any([k,n] == 0)
		return;
	end
	
	iunkvars = setdiff((1:k)', iknownvars);
	Y(iunkvars, end) = 0;
	Z(iunkvars, :) = 0;
	
	J = Op.Jmat;
	JR = zeros(k);
	JR(:, iknownvars) = (-J(:, iknownvars));
	JR = JR - diag(diag(JR));
	JL = J + JR;
	
	p = Op.Order();
	if first > p
		Y = Y(:, first-p:last);
	else
		padlen = p-first+1;
		Y = [zeros(k, padlen), Y(:, 1:last)];
	end
	
	ll = Op.lags;
	Aunk = Op.Phi(iunkvars, :, :);
	if Op.has_lag_zero
		ll = [0, ll];
		Aunk = cat(3, JR(iunkvars, :), Aunk);
	end
	
	% Now apply operator to Y.
	% Code vectorized over index range, first to last.
	Y(~isfinite(Y)) = 0; % Want 0*NaN and 0*Inf = 0.
	for i = 1 : numel(ll)
		Z(iunkvars, :) = Z(iunkvars, :) + Aunk(:, :, i)*Y(:, p+1-ll(i) : end-ll(i));
	end
	
	Z(iunkvars, :) = JL(iunkvars, iunkvars)\Z(iunkvars, :);
end
%%
