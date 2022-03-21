% ParseArgs
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
function [Lags, Coeffs, J, has_lag_zero] = ParseArgs(Lags, Coeffs, J)
	ERR_COEFFSIZES = ['Lag Operator:Coefficient array must have a shape that ', ...
		'is consistent with an array of square matrix slices. '];
	
	k = sqrt(numel(Coeffs)/numel(Lags));
	assert(rem(k,1) == 0, ERR_COEFFSIZES);
	
	sparse = issparse(Coeffs);
	if sparse
		Coeffs = full(Coeffs);
	end
	
	% Bring coefficients into 3D format.
	Coeffs = reshape(Coeffs, k, k, numel(Lags));
	
	assert(isnumeric(Coeffs) && (ndims(Coeffs) <= 3), ERR_COEFFSIZES);
	
	[k, l, q] = size(Coeffs);
	assert(k > 0, 'LagOperator:To construct an empty operator, use Either a (k x k x 0) array Or specify all-zero coefficients. ');
	assert(l == k, ERR_COEFFSIZES);
	assert(isempty(Lags) || isvector(Lags));
	assert(all(Lags > 0));
	assert(numel(Lags) == q, 'Lag Operator:Number of lags does not match number of coefficients.');
	Lags = vecr(Lags);
	
	if isempty(J)
		J = eye(k);
		has_lag_zero = false;
		
	else
		assert(isnumeric(J) && ismatrix(J));
		[k1, k2] = size(J);
		assert((k1 == k) && (k2 == k));
		assert(all(diag(J) == 1), 'All diagonal elements of J must be equal to 1.');
		has_lag_zero = ~isequal(J, eye(k));
	end
	
	% Sort lags and coefficients.
	if ~issorted(Lags)
		[Lags, ixvec] = sort(Lags);
		Coeffs = Coeffs(:, :, ixvec);
	end
	
	% Remove zero coefficients.
	[Lags, Coeffs] = Trim(Coeffs, Lags);
	
	if sparse
		Coeffs = sparse(Coeffs);
	end
end
%%
