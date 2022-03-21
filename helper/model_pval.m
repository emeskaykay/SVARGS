% model_pval
% Returns p-value(s) for the hypothesis that the log-likelihood ratio is greater than 0.
%
% For general error covariance matrices analogs of the F-statistic are the likelihood ratio
% or the granger causality. The likelihood ratio reduces to (1 + const*F_statistic) for the
% case of scalar error covariance. Unfortunately, the distribution of the likelihood ratio
% in the general case is not commonly known. The workaround is to take -2*log_likelihood_ratio
% which, by Wilk's theorem, has, asymptotically, a known Chi-square distribution.
%
% detsig_0:
% One of:
%	Scalar error variance of the original model.
% OR
%	Vector of scalar error variances of the original model.
% OR
%	Error covariance matrix of the original model.
%
% detsig:
%	Scalar error variance of the updated model.
% OR
%	Vector of scalar error variances of the updated model.
% OR
%	Error covariance matrix of the updated model.
%
% cxdelta:
%	Complexity of updated model - Complexity of original model.
%	If the 'direction' input is provided, the sign of cxdelta is ignored.
%
% datalen:
%	The number of data points that were used to compute the maximum likelihoods (the maximum
%	likelihoods are proportional to 1./detsig_0 and 1./detsig).
%
% direction:
%	signed value.
%	+1 indicates that the complexity of the updated model is greater than that of the original.
%	-1 indicates that the complexity of the updated model is less than that of the original.
%	0 indicates that the model complexities cannot be compared. The corresponding pvalues will
%	be NaN.
%
% NOTE:
%	detsig_0, detsig and cxdelta can all be vectors of the same length and orientation, in
%	which case pvals is the corresponding vector of p-values. The arguments datalen and
%	direction can, each, either be scalars or vectors of the same length as the first three
%	arguments.
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
function [pvals, gcstrength] = model_pval(detsig_0, detsig, cxdelta, datalen, direction)

	if (nargin < 5) || isempty(direction)
		direction = 2*( (cxdelta > 0) - 1/2 );
	else
		cxdelta = abs(cxdelta);
	end
	direction = sign(direction);
	
	if ~isvector(detsig_0)
		assert( is_square(detsig_0) && all(size(detsig_0) == size(detsig)) );
		sig_0 = detsig_0;
		sig = detsig;
		
		if rank(sig_0) < length(sig_0)
			gcstrength = 0;
			
		else
			logdetsig_0 = log_determinant(sig_0); % log of absolute value of determinant
			logdetsig = log_determinant(sig);
			gcstrength = direction*(logdetsig_0 - logdetsig);
		end
		
	else
		if isscalar(direction), direction = direction*ones(size(detsig_0)); end
		if isscalar(datalen), datalen = datalen*ones(size(detsig_0)); end
		
		gcstrength = zeros(size(detsig));
		
		bnz = (detsig_0 > 0);
		gcstrength(bnz) = direction(bnz).*( log(detsig_0(bnz)) - log(detsig(bnz)) );
	end
	
	llr = 0.5*datalen.*gcstrength;
	
	% chi2cdf(x, v, 'upper') returns the complement of the chi-square cdf at each value in x,
	% using an algorithm that more accurately computes the extreme upper tail probabilities.
	pvals = chi2cdf(2*llr, cxdelta, 'upper');
	assert(all(pvals >= 0) && ~any(isinf(pvals)));
	
end
%%
