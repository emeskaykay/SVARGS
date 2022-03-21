% [Op, Eigs] = var_rand_operator(k, p, rlims)
% Generates a 'random' stable VAR operator, given the number of variables
% and the order of the operator. The operator is generated by choosing
% eigenvalues and eigenvectors randomly. See the functions:
%	rand_disk()
%	var_rand_coeffs()
% The number of complex eigenvalue pairs is chosen to be a random number
% between 0 and k*p/2. The remaining eigenvalues are real.
% 
% k: The number of variables that the operator must contain.
%
% p: The order of the operator.
%
% rlims [default [0, 1]):
%	The limits on the magnitude of the eigenvalues. This can also be a matrix
%	of size ((k*p - q) x 2) giving individual limits for each of the k*p eigenvalues.
%	Here q is the number of complex eigenvalue pairs. Only one of each the complex
%	conjugate pairs can have their limits specified.
%
% Returns:
% Op:
%	A LagOperator object. See class LagOperator.
%
% Eigs:
%	The eigenvalues of the model.
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
function [Op, Eigs] = var_rand_operator(k, p, rlims)

	% Total number of Eigenvalues.
	nEigs = k*p;
	if (nargin < 3) || isempty(rlims), rlims = [0,1]; end
	assert(size(rlims, 2) == 2);
	
	if numel(rlims) == 2
	% Number of pairs of complex e-values.
		q = floor(0.5*nEigs*rand);
	else
		q = k*p - size(rlims,1);
		assert(q >= 0);
	end

	% Generate (nEigs-q) random Eigenvalues of which q are uniformly
	% distributed in the open unit disk in the (x,y) plane, and the
	% remaining nEigs-2*q are real and uniformly distributed on [0,1].
	ttalims = [[-pi*ones(q, 1), pi*ones(q, 1)]; zeros(nEigs-2*q, 2)];
	Eigs = rand_disk([zeros(nEigs-q,1), ones(nEigs-q,1)], ttalims);

	% Adjust the eigenvalue lengths to get a bell shaped, (rather than
	% uniform) distribution for the square of the length. We use a
	% beta(3,3) distribution here.
	r = abs(Eigs);
	rsqr = r.^2;
	rsqradj = icdf('beta', rsqr, 3, 3);
	Eigs = sqrt(rsqradj./rsqr).*Eigs;
	r = abs(Eigs);

	% Squeeze the eigenvalues into the specified ranges.
	Eigs = ( rlims(:,1) + (rlims(:,2) - rlims(:,1)).*r ).*(Eigs./r);

	% Add the complex conjugates of the non real eigenvalues.
	conj_eigs = conj(Eigs(Eigs ~= real(Eigs)));
	Eigs = sort([Eigs; conj_eigs]);

	% Generate Coefficients.
	C = var_rand_coeffs(k,p,Eigs);
	% ___

	Op = LagOperator(1:p, C);
end
%%