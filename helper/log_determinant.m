% Computes the natural log of the absolute value of the determinant
% of a matrix without computing the actual determinant. This is useful
% when there is significant overflow or underflow when computing the
% determinant itself.
%
% A
%	Square matrix.
%
% spd [default: true]:
%	When spd is true, the matrix is (initially) assumed to be symmetric
%	positive definite.
%
% RETURNS:
% d:
%	Effectively log(|det(A)|).
%
% spd:
%	true if A is actually symmetric positive definite, false otherwise.
%
% NOTES:
%	Uses Cholesky or QR factorization.
%
% NOTE 2:
%	The result is the log of the absolute value of the determinant, so
%	it is always real even if one or both of the matrices are complex.
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
function [d, spd] = log_determinant(A, spd)

	if (nargin < 2) || isempty(spd), spd = true; end
	
	if isscalar(A)
		detA = A;
		spd = (isreal(detA) && (detA > 0));
		d = log(abs(detA));
		
	elseif spd
		% note that spd => all diagonal elements of R will be real positive.
		[R,p] = chol(A);
		if p == 0
			r = diag(R);
			d = 2*sum(log(r));
		else
			spd = false;
			d = log_determinant(A, spd);
		end
	
	else
		% triu(X) is the upper-triangular matrix, R, of the QR factorization of A.
		X = qr(A);
		r = abs(diag(X));
		d = sum(log(r));
	end

end
%%
