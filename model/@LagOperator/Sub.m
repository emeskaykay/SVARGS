% Sub
% returns Sub-operator.
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
function [Op, kseq] = Sub(Op, kseq, lagseq)
	
	if (nargin < 3), lagseq = ':'; end
	
	if islogical(kseq), kseq = find(kseq); end
	if islogical(lagseq), lagseq = find(lagseq); end
	
	k = NumVars(Op);
	
	% If kseq is the string 'nontriv', then elimnate those components
	% that have no dependencies with other components. (ie: their
	% order is zero and no other component depends on them.
	if ischar(kseq) && strncmpi(kseq, 'nontriv', 7)
		C = Coeffs(Op, Lags(Op));
		ktriv_col = all(all(C == 0, 2), 3);
		ktriv_row = all(all(C == 0, 1), 3);
		ktriv = (ktriv_col & ktriv_row');
		if ~isempty(ktriv)
			krange = (1 : k);
			kseq = krange(~ktriv);
		else
			kseq = ':';
		end
	end
	
	assert( isvector(lagseq) && isvector(kseq) );
	
	allk = is_colon(kseq) ...
		|| ( (length(kseq) == k) && all(kseq == 1:k) );
	alllags = is_colon(lagseq) ...
		|| ( (length(lagseq) == length(Op.lags)) && all(lagseq == Op.lags) );
	
	if (allk && alllags), return; end % Original operator returned.
	% ----------- %
	
	if (allk), kseq = 1:k;
	elseif (alllags), lagseq = Op.lags;
	end
	
	assert( is_unique(kseq) );
	assert( is_unique(lagseq) );
	
	assert( length(kseq) <= NumVars(Op) );
	assert( length(lagseq) <= Order(Op) );
	
	[C, ix, J] = Coeffs(Op, lagseq);
	
	Csub = C(kseq,kseq,ix);
	Jsub = J(kseq,kseq);
	lagssub = lagseq(ix);
	[Csub, lagssub] = Trim(Csub, lagssub);
	
	Op = LagOperator(lagssub, Csub, 'J', Jsub);
	
end
%%
