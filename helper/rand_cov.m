% rand_cov
% Generates a random covariance matrix.
%
% k:
%	The size of the covariance matrix.
%
% dx [default: 0.1 + rand]
%	Either a positive scalar specifying the determinant of the covariance matrix,
%	OR
%	A positive vector of length k, specifying the diagonal (variances) of the covariance matrix.
%
% varargin:
% The covariance matrix is generated via a random unitary matrix.
% Name value pair options for generating random unitary matrix can be
% specified as below:
%	'uniform_cosines', value: true/false [default: true]:
%		If 'uniform_cosines' is true (faster), the direction cosines of the axes
%		(rather than the angles) are uniformly distributed. If false, it is the
%		angles that are uniformly distributed.
%
% Returns:
% Sigma:
%	A random k x k covariance matrix (ie: Symmetric positive-definite matrix)
%	with the determinant OR diagonal as specified in the dx argument.
%
% detval:
%	Determinant of Sigma.
%
% ll:
%	Eigenvalues of Sigma.
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
function [Sigma, detval, ll] = rand_cov(k, dx, varargin)

	if (nargin < 2) || isempty(dx), dx = 0.1 + rand; end
	assert(isnumeric(dx) && isvector(dx));
	
	if isscalar(dx)
		assert(dx > 0);
		detval = dx;
		diagvals = [];
		
	else
		assert(numel(dx) == k);
		assert(all(dx > 0));
		diagvals = dx(:);
		detval = 1;
	end
	
	if k == 1
		Sigma = dx;
		detval = dx;
		ll = dx;
		
	else
		% Get random unitary matrix.
		U = rand_unitary(k, varargin{:});

		% Generate covariance matrix.
		ll = generate_evalues__(k, detval);
		Sigma = U*diag(ll)*U';

		if ~isempty(diagvals)
			% Normalize Sigma so that the ratio of the diagonal entries are equal to 'diagvals'.
			Sigma = Normalize( Sigma, diag(diag(Sigma)./diagvals) );
		end

		Sigma = (Sigma + Sigma')/2;
	end
	
end
%%

function ll = generate_evalues__(k, prodval)
	ll = 0.1 + rand(k,1);
 	xinvert = datasample(1:k, floor(k/2));
 	ll(xinvert) = 1./ll(xinvert);
	logr = (1/k)*( log(prodval) - sum(log(ll)) );
	ll = ll*exp(logr);
end
%%
% 
