% Evalues
% Usage: Evalues(Resolution)
% Returns the eigenvalues of the VAR(p) operator:-
%		A0.y(t) - A1.y(t-1) - A2.y(t-2) - ... - Ap.y(t-p)
%	=	(A0 - Phi(L))(Y(t))
%
% Resolution [Optional. See Eigenvalues.m for defaults]:
%	Tolerance for eigenvalues to be treated as equal (See Eigenvalues.m).
%
% Returns:
% l:
%	Sorted column vector of eigenvalues of the characteristic polynomial.
%
% ldistinct:
%	Sorted column vector of distinct eigenvalues of the characteristic polynomial.
%
% lmultiplicity:
%	Vector of multiplicities corresponding to the eigenvalues in ldistinct.
%
% Note that the eigenvalues returned are those of the charcteristic polynomial
% of the associated AR(1) process, not the 'reverse' characteristic polynomial.
% To get the eigenvalues of the reverse characteristic, take the reciprocal.
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
function [l, ldistinct, lmulti] = Evalues(Op, Resolution)

	k = NumVars(Op);
	p = Order(Op);
	if isempty(Op.lags)
		l = 0;
		
	else
		if nargin < 2, Resolution = []; end
		
		dataid = ['evalues_', num2str(Resolution)];
		[l, ldistinct, lmulti] = GetSavedEvalues__(Op, dataid);
		if ~isempty(l)
			return;
		end
	
		if ~MemoryCheckMsg(2*(k*p)^2)
			error('LagOperator::Evalues: Execution cancelled. ');
		end

		A = BlockCompanion(Op.Coeffs(1:p, 'reduced', true));
		if (nargout <= 1)
			l = Eigenvalues(A, [], 'Resolution', Resolution);
		else
			[l, ldistinct, lmulti] = Eigenvalues(A, [], 'Resolution', Resolution);
		end
		
		SaveEvalues__(Op, dataid, l, ldistinct, lmulti)
	end
	
end
%%

%%
function [l, ldistinct, lmulti] = GetSavedEvalues__(Op, dataid)
	data = Op.hCache.Get(dataid);
	l = [];
	ldistinct = [];
	lmulti = [];
	if ~isempty(data)
		l = data.l;
		ldistinct = data.ldistinct;
		lmulti = data.lmulti;
	end
end

%%
function SaveEvalues__(Op, dataid, l, ldistinct, lmulti)
	data = Op.hCache.Get(dataid);
	data.l = l;
	data.ldistinct = ldistinct;
	data.lmulti = lmulti;
	Op.hCache.Set(dataid, data);
end
%%
