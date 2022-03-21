% val = sbc_x(mosobj, cost, m)
% Extended Schwartz Bayesian Criterion value.
% For cases when M is comparable to N or or larger than N.
%
% cost:
%	log absolute determinant of error covariance. This can be a vector of values.
%
% m:
%	Total number of parameters (complexity).
%	This can either be a scalar or a vector of values of the same length as 'cost'.
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
function val = sbc_x(mosobj, cost, m)

	N = mosobj.N;
	M = mosobj.M;
	k = mosobj.k;
	
	sbcval = sbc(mosobj, cost, m);
	
	T = zeros(size(m));
	for i = 1 : numel(m)
		r = 0 : m(i)-1;
		T(i) = sum(log(M - r) - log(r + 1));
	end
	
	% M == N^kappa. Want:
	% gamma ~= 0.5 for N well above M
	% gamma ~= 1 for N well below M
	z = inv_std_normal(1 - 0.01/(2*k*M));
	kappa = log(M)/log(N);
	a = normcdf(-log2(kappa), 0, 1/z);
	gamma = 1 - a/2;
	val = sbcval + 2*gamma*T/N;
end
%%
