% qr_factor
% Wrapper around Matlabs qr(X,0). In addition the function returns the indices of
% the regressor columns, if any, that are linearly dependent.
%
% X:
%	An (N x m) matrix.
%
% pivot_tol [default: 0]:
%	If pivot_tol > 0, the rank revealing QR factorization is used which implies
%	column pivoting. This will result in a slower, but possibly more stable
%	factorization. In addition, the outputs, xperm and depcount, may have non-trivial
%	values - see below.
%
% Returns:
% R:
%	An (min(m,N) x m) upper triangular matrix.
%
% Q:
%	An (N x min(m,N)) unitary matrix.
%
% xperm (Required only if pivot_tol > 0):
%	This is a (1 x m) permutation vector such that:
%			X(:, xperm)   =   X*E   =   Q*R
%		where E = I(:,xperm)
%
%		which is the same as:
%			X   =   Q*R(:, invperm)   =   Q*R*E'
%		(where invperm is the inverse permutation of the cperm vector)
%
% depcount (Required only if pivot_tol > 0):
%	This is the number of dependent columns of X (m - rank(X)).
%	The dependent columns of X are unused. The indices of the dependent
%	columns of X are placed in the last depcount elements of xperm:
%		ie: ixdep = xperm(end-dpecount+1 : end);
%
% If pivot_tol == 0, then the last two outputs are (1:size(R,2)) and 0 respectively.
%
% pivot_tol > 0 will result in a slower, but possibly more numerically stable
% factorization using column pivoting.
%
% Specifying only one output (R) will lead to the fastest return.
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
function [R, Q, xperm, depcount] = qr_factor(X, pivot_tol)

	if (nargin < 2) || isempty(pivot_tol), pivot_tol = 0; end
	
	if pivot_tol == 0
		if (nargout <= 1)
			R = qr(X, 0);
			
		else
			[Q, R] = qr(X, 0);
			xperm = (1 : size(R,2));
			depcount = 0;
		end
		
	else
		assert(nargout >= 3, 'Need 3 output args (QR with pivoting returns column permutation). ');
		
		% Use the rank-revealing QR to indicate dependent columns of X.
		[Q, R, xperm] = qr(X, 0);
		
		[N,m] = size(X);
		Xp = X(:, xperm(1:min(m,N)));
		ss = vecc(sqrt(sum(Xp.*Xp, 1)));
		rdiag = abs(diag(R));
		l = numel(rdiag);
		l0 = find(rdiag >= pivot_tol*ss, 1, 'last');
		depcount = max(0, m-N) + l - l0;
	end

end
%%
