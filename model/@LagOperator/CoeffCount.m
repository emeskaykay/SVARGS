% CoeffCount
% Returns number of non-zero coefficients.
%
% varargin:
% Name-Value Pairs;
%	Name: 'reduced', Value: true/false. [default: false]
%	If reduced is true, the coefficient counts are that of the reduced model.
%	(See LagOperator::Reduce.m)
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
function [n, ijcounts] = CoeffCount(Op, ixeq, varargin)

	if nargin < 2, ixeq = ':'; end
	
	k = Op.NumVars();
	bnr = true;
	if Op.has_lag_zero
		opts = parse_nvpairs(varargin, {'reduced', 'logical', false});
		bnr = (false == opts.reduced);
	end
	
	if nargout == 1
		if is_colon(ixeq) || isequal(ixeq(:), (1:k)')
			n =  nnz(Op.Phi___) + bnr*(nnz(Op.Jmat__) - k);
			
		else
			q = numel(Op.lags);
			ix = vecc(ixeq) * (1:q);
			n =  nnz(Op.Phi___(ixeq, vecr(ix))) + bnr*(nnz(Op.Jmat__(ixeq, ixeq)) - numel(ixeq));
		end
		
	elseif nargout > 1
		dataid = {'eq_coeff_counts', bnr};
		eqdata = Op.hCache.Get(dataid);
		if isempty(eqdata)
			B = logical(Op.Phi___);
			B0 = bnr * (logical(Op.Jmat__) & ~logical(eye(k)));
			eqdata.counts = full(sum(B(:,:), 2) + sum(B0,2));
			eqdata.total_count = sum(eqdata.counts);
			Op.hCache.Set(dataid, eqdata);
		end
		ijcounts = eqdata.counts;
		n = sum(ijcounts(ixeq));
	end
	
end
%%
