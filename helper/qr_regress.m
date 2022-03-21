% qr_regress
% Regression based on QR factorization of regressor matrix. This is useful when
% there are multiple (k) response variables for the same set of regressors. We
% can then do just a single QR factorization, which would otherwise have to be
% performed k times. This speeds up the computation immensely.
%
% Y:
%	A (N x k) matrix where N is the length if the time series. The i-th column
%	is the response vector for the i-th component.
%
% X:
%	(N x m) matrix of regressors. The regressors must be the columns of X.
%
% breglocs [ default: [] ]:
%	Boolean row vector mask of length m.
%	breglocs indicates the regressors (columns of X) that are to be used to fit
%	the response variable Y. If breglocs is empty, all columns of X are used.
%
% include_intercept [optional] [default: false]:
%	If true, a constant ones regressor is added to the set of regressors. The
%	output then contains a non-zero intercept row vector (c0). If false, then c0
%	is the zero row vector. Note that in either case there no check to see whether
%	X already has a column of ones.
%
% pivot_tol [optional] [default: 0]:
%	If pivot_tol > 0, the rank revealing QR factorization is used which implies
%	column pivoting. This results in a slower, but possibly more stable regression.
%	In addition, the output ixdep may be non-empty - see below.
%
% Returns:
% C:
%	An (m x k) coefficient matrix where k = size(Y,2) is the number of components and
%	m = size(X,2) is the number of regressors. X*C is the least squares fit for Y and
%	Res = Y - X*C are the residuals. If ixdep is non-empty (which is possible only if
%	pivot_tol > 0), then C has the rows correponding to ixdep set to 0.
%
% c0:
%	If include_intercept is true, a ones column is appended to the set of regressors
%	and c0(i) is the coefficient of the ones regressor (intercept) for the i-th component
%	of Y. If include_intercept is false, c0 is the zero row vector.
%
% R:
%	An (m x m) upper triangular matrix in the QR factorization of X.
%	If ixdep is non-empty (which is possible only if pivot_tol > 0), then R(ix,ix) is
%	the R term in the QR factorization of X(:, ix), where ix = setdiff(1:m, ixdep), and
%	the rows and columns of R corresponding to ixdep are set to 0. R'*R is equal to the
%	covariance X'*X for all pairs of indices (i,j) for which i and j are each in ix.
%	For (i,j) in ixdep, R'*R is identically 0.
%
%	Note that this is for the factorization _without_ the intercept regressor. Since
%	the column of ones, if added, is appended as the last column of X, this just
%	involves removing the last column of Q and the last row and column of R from the
%	QR factorization of [X, ones].
%
% ixdep (only if pivot_col > 0):
%	Indices of the "dependent" (upto pivot_tol) columns of X. If pivot_tol == 0, then
%	ixdep is empty.
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
function [C, c0, R, ixdep] = qr_regress(Y, X, breglocs, include_intercept, pivot_tol)
	[N, m] = size(X);
	[~, k] = size(Y);
	assert(size(Y,1) == N);
	
	if (nargin < 5) || isempty(pivot_tol), pivot_tol = 0; end
	if (nargin < 4) || isempty(include_intercept), include_intercept = false; end
	if (nargin < 3) || isempty(breglocs), breglocs = true(1,m); end
	
	assert(isscalar(include_intercept) && islogical(include_intercept));
	assert( isvector(breglocs) && islogical(breglocs) && (numel(breglocs) == size(X,2)) );
	X = X(:, breglocs);
	clear breglocs
	
	if include_intercept
		X = [X, ones(N,1)];
		m = m+1;
	end
	
	c0 = zeros(1,k);
	ixdep = [];
	if N >= m
		[C, R, ~, ixdep] = lsfit__(X, Y, pivot_tol);
		
	else
		% Option 1: Take Pseudo-inverse.
		% Option 2: Keep only the N most orthogonal columns of X
		% Currently we have chosen option 1.
		C = pseudoinv__(X, Y, pivot_tol);
		if nargout > 2
			[R, ~, xperm, depcount] = qr_factor(X, pivot_tol);
			ixdep = xperm(end-depcount+1:end);
			
			% Do this to keep most orthogonal columns of X
% 			C = zeros(m,k);
% 			ix = setdiff(1:m, ixdep);
% 			[C(ix, :), R] = lsfit__(X(:, ix), Y, 0); % Unpermuted
		end
	end
	
	if include_intercept
		% Save the coefficients of the ones column in c0 and remove the same from C.
		c0 = C(end, :);
		C(end, :) = [];

		if nargout > 2
			R(end,:) = [];
			R(:,end) = [];
		end
	end
end

%%
% QR factorization of X:
% Q*R = X(:, xperm)
% Q is (N x m)
% R is (m x m)
% X is (N x m)
%
% For ls solution:
%		X*C = Y:
% =>	X(:, xperm)*C(xperm, :) = Y
% =>	X(:, xperm)'*X(:, xperm)*C(xperm, :) = X(:, xperm)'*Y
% =>	R'*Q'*Q*R*C = R'*Q'*Y
% =>	R'*R*C = R'*Q'*Y
% =>	R*C = Q'*Y
% =>	C = R \ Q'*Y
% Now if we omit the last depcount elements of xperm, then letting
%		m0 = m - depcount
%		xperm0 = xperm(1:m0)
%		R0 = R(:, 1:m0)
% the above calculation becomes:
%		C(xperm0, :) = = R0 \ (Q'*Y);
%
% If there is no pivoting we can also use the following lines of code:
%		R = qr_factor(X, 0);
% 		Q = [];
% 		R = triu(R(1:m, :));
%		C = R \ (R' \ (X'*Y))
% which has the advantage of not needing the matrix Q. This helps
% since matlabs S = qr(A,0) is about twice as fast as [Q,R] = qr(A,0).
% However the two inversions cause the resulting solution to be less
% accurate.
function [C, R, Q, ixdep] = lsfit__(X, Y, pivot_tol)

	if pivot_tol == 0
		% No pivoting
		[R, Q] = qr_factor(X, 0);
		C = R \ (Q'*Y);
		ixdep = [];
		
	else
		[R, Q, xperm, depcount] = qr_factor(X, pivot_tol);
		
		[~, k] = size(Y);
		[~, m] = size(X);
		C = zeros(m, k);
		ixdep = [];
		if depcount > 0
			ixdep = xperm(end-depcount+1:end);
			ix = setdiff(1:m, ixdep);
			[R(ix, ix), Q(:, ix)] = qr_factor(X(:, ix), 0); % Unpermuted
			R(ixdep, :) = 0;
			R(:, ixdep) = 0;
			Q(:, ixdep) = 0;
			
			% Compute the LS coefficients.
			C(ix, :) = R(ix,ix) \ (Q(:,ix)'*Y);

		else
			% Compute the LS coefficients.
			C(xperm, :) = R \ (Q'*Y);
			
			% No pivoting
			R = qr_factor(X, 0);
		end
	end
	
end

%%
% N < m. Underdetermined system. Moore-Penrose pseudoinverse.
% This gives the minimum norm solution for each column of C.
% Use QR factorization of transpose(X):
%		X(xperm, :)' = Q*R
%
% Then to solve:
%		X*C = Y
%
% use the Moore-Penrose preudoinverse. By definition:
%		C = X'*inv(X*X')*Y
%
% =>	C = X(xperm, :)'*inv(X(xperm, :)*X(xperm,:)') * Y(xperm, :)
% =>	C = Q*R*inv((Q*R)'*(Q*R))*Y(xperm, :)
%		C = Q*R*inv(R'*Q'*Q*R)*Y(xperm, :)
% =>	C = Q*R*inv(R'*R)*Y(xperm, :)
% =>	C = Q * (R' \ Y(xperm, :))
%		
% If there is no pivoting we can also use the following two lines of code:
%		R = qr_factor(Xt, 0);
%		C = Xt * (R \ (R' \ Y));
% which has the advantage of not needing the matrix Q. This helps
% since matlabs S = qr(A,0) is about twice as fast as [Q,R] = qr(A,0).
% However the two inversions cause the resulting solution to be less
% accurate.
function C = pseudoinv__(X, Y, pivot_tol)

	Xt = transpose(X);
	if pivot_tol == 0
		% No pivoting
		[R, Q] = qr_factor(Xt, 0);
		C = Q * (R' \ Y);
		
	else
		[R, Q, xperm, depcount] = qr_factor(Xt, pivot_tol);
		
		if depcount == 0
			C = Q * (R' \ Y(xperm, :));
			
		else
			[~, m] = size(X);
			[~, k] = size(Y);
			C(m, k) = 0; % Implicit Initialization.
			ixdep = xperm(end-depcount+1:end);
			ix = setdiff(1:m, ixdep);
			[R(ix, ix), Q(:, ix)] = qr_factor(Xt(:, ix), 0); % Unpermuted
			C(ix, :) = Q(:,ix) * (R(ix,ix)' \ Y);
		end
	end
	
end
%%
