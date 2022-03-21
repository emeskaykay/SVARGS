% Compose__
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
function Op = Compose__(leftOp, rightOp)
	% If either operator is empty, return the other operator
	if isempty(leftOp.lags) && (false == leftOp.has_lag_zero)
		Op = rightOp;
		return;
		
	elseif isempty(rightOp.lags) && (false == rightOp.has_lag_zero)
		Op = leftOp;
		return;
	end
	
	k = NumVars(leftOp);
	I = eye(k);
	leftlags = leftOp.lags;
	rightlags = rightOp.lags;
	ql = length(leftlags);
	qr = length(rightlags);
	[Cleft, ~, Jleft] = leftOp.Coeffs(leftlags);
	[Cright, ~, Jright] = rightOp.Coeffs(rightlags);
	
	if leftOp.has_lag_zero
		Cleft = cat(3, Jleft - I, Cleft);
		leftlags = [0, leftlags];
		ql = ql + 1;
	end
	
	if rightOp.has_lag_zero
		Cright = cat(3, Jright - I, Cright);
		rightlags = [0, rightlags];
		qr = qr + 1;
	end
	
	% Calculate all possible pairwise sums of the lags.
	Lagmat = zeros(ql, qr);
	for col = 1:qr
		for row = 1:ql
			Lagmat(row,col) = leftlags(row) + rightlags(col);
		end
	end
	
	% Full column linear sort of matrix. 
	[Lags, linix] = sort(Lagmat(:));

	% Get the subscripted (i,j) indices of the linear indices.
	[ixvec, jxvec] = ind2sub([ql qr], linix);
	
	qres = leftOp.has_lag_zero + rightOp.has_lag_zero + leftlags(end) + rightlags(end);
	qres = min(ql*qr, qres);
	UniqueLags = zeros(1, qres);
	Cres = zeros(k, k, qres);
	
	n = 1;
	UniqueLags(1) = Lags(1);
	for c = 1 : ql*qr
		% Note that Lags is already sorted so that Lags(c) cannot be less than Lags(n).
		if (Lags(c) > UniqueLags(n))
			n = n+1;
			UniqueLags(n) = Lags(c);
		end
		
		% Convert (row, col) to the corresponding (i,j) which ensures that the 
		% lag sums corresponding to Phi(:,:,i)*Phi(:,:,j) come up in sorted order.
		i = ixvec(c);	j = jxvec(c);
		
		Cres(:,:,n) = Cres(:,:,n) + Cleft(:,:,i)*Cright(:,:,j);
	end
	
	% Blank out unused slots in Lags and Phi.
	UniqueLags(n+1:end) = [];
	Cres(:,:,n+1:end) = [];
	
	if UniqueLags(1) > 0
		Op = LagOperator(UniqueLags, Cres);
		
	else
		J = I + Cres(:,:,1);
		Cres = Cres(:,:,2:end);
		UniqueLags = UniqueLags(2:end);
		Op = LagOperator(UniqueLags, Cres, 'J', J);
	end
end
%%
