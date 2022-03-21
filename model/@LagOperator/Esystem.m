% Esystem
% Usage: Esystem(Params)
% Returns the Eigenvalues and eigenvectors of the characteristic polynomial
% of the VAR(p) operator:-
%			A0.y(t) - A1.y(t-1) - A2.y(t-2) - ... - Ap.y(t-p)
%		=	( A0 - Phi(L) )(Y(t))
%
% Note that the eigenvalues returned are those of the charcteristic polynomial
% of the associated AR(1) process, not the 'reverse' characteristic polynomial.
% To get the eigenvalues of the reverse characteristic polynomial, take the reciprocal.
%
% Params:
%	Eigenparams object containing eigensystem parameters (See Eigenparams.m)
%
% Return values:
% Q:
%	A kxk matrix with the eigenvectors forming the columns of Q. The left to
%	right ordering of eigenvectors in Q correspond to the diagonal ordering in L.
%
% L:
%	A kxk diagonal matrix of eigenvalues.
%
% Code:
%	See comments in Eigensystem.m
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
function [Q L Code] = Esystem(Op, Params)
	
	if (isempty(Op.lags))
		L = []; Q = []; Code = []; return; end
	
	if (nargin < 2) || isempty(Params), Params = Eigenparams(); end
	assert(isa(Params, 'Eigenparams'));
	
	dataid = {'esystem', Params};
	[Q L Code] = GetSavedEsystem__(Op, dataid);
	if ~isempty(Q)
		return;
	end
	
	k = NumVars(Op);
	p = Order(Op);
	n = k*p;
	
	% Note: Evalues returns sorted eigenvalues.
	[evec, edistinct, emulti]= Evalues(Op, Params.Resolution);
	assert(length(evec) == n);
	
	numdistinct = length(edistinct);
	if numdistinct == n
		[Q, L] = eig(BlockCompanion(Op.Coeffs(1:p)));
		Code = 'DIAGONAL';
		SaveEsystem__(Op, dataid, Q, L, Code);
		return;
	end
	
	disp('Esystem:Operator has ''repeated'' eigenvalues. ');
	
	nsdim = ones(numdistinct, 1); % vector of nullspace dimensions
	I = eye(k);
	Iarray = repmat(I, [1,1,p-1]);
	last = 0;
	lags = Lags(Op);
	Q = zeros(n,n); % Eigenvector matrix.
	C = Coeffs(Op, lags, 'reduced', true);
	for i = 1 : numdistinct
		e = edistinct(i);
		lambda_array = cumprod(e*Iarray, 3);
		lambda_array = lambda_array(:, :, end:-1:1);
		lambda_array = cat(3, lambda_array, I);
		
		% Form the matrix E.
		% E = e^p*I - e^(p-1)*A(1) - e^(p-2)*A(2) - ... - e^1*A(p-1) - A(p).
		E = e^p - DotM(C, lambda_array(:, :, lags));
		
		% Find the null space of E.
		m = emulti(i);
		Z = NullSpace(E, Params.rankTol, 1, min(m, k));
		
		% Form the matrix [e^(p-1)*Z; e^(p-2)*Z; . . . ; e^1*Z; Z]
		% and assign it to the correct columns of Q.
		nsdim(i) = size(Z, 2);
		lambda_mat = reshape(lambda_array, n, k);
		Q(:, last+1 : last+nsdim(i)) = lambda_mat*Z;
		
		last = last + m;
	end
	
	if ( sum(nsdim) == n )
		Code = 'DIAGONAL';
		
	else
		% At this point operator is defective.
		% Construct the full AR(1) matrix.
		A = BlockCompanion(Op.Coeffs(1:p));
		
		Code = 'DEFECTIVE';
		disp('Eigensystem: Matrix is defective. ');
		
		if Params.bExactify
			disp('Eigensystem: Exactifying defect. ');
			A = ExactifyDefect(A, evec, Q);
		end

		if (Params.bExactify && Params.bJordan)
			% At this point the matrix A (after exactification) is exactly defective.
			disp('Eigensystem: Computing Jordan form. ');
			disp('Eigensystem: Assuming numerical values are ''very'' accurate. ');
			try
				Q = jord(A, evec, Params.rankTol); % ExtSrc__
				Code = 'JORDAN';
				
			catch err
				disp('Eigensystem: Jordan Form Computation Failed. ');
				disp(err.message);
				Code = 'FAILED';
			end
		end
	end
	L = diag(evec);
	
	SaveEsystem__(Op, dataid, Q, L, Code);
end

%%
function [Q L Code] = GetSavedEsystem__(Op, dataid)
	data = Op.hCache.Get(dataid);
	Q = [];
	L = [];
	Code = [];
	if ~isempty(data)
		Q = data.Q;
		L = data.L;
		Code = data.Code;
	end
end

%%
function SaveEsystem__(Op, dataid, Q, L, Code)
	data = Op.hCache.Get(dataid);
	data.Q = Q;
	data.L = L;
	data.Code = Code;
	Op.hCache.Set(dataid, data);
end
%%
