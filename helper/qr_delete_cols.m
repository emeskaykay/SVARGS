% qr_delete_cols
% Can use full (m x m) or 'economy' sized Q (m x n)
% R must correspondingly be full (m x n) or economy sized (n x n).
%
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
function [R, Q] = qr_delete_cols(ixcols, R, Q)

	assert(nargin == 3, 'qr_delete_cols: 3 arguments required. ');
	
	m = size(Q,1);
	n0 = size(R,2);
	p = numel(ixcols);
	
	assert(size(Q,2) == size(R,1));
	assert(all(ixcols <= n0));
	
	if p == n0
		R = zeros(0,0);
		Q = zeros(m,0);
		return;
	end
	
	% More efficient if higher index columns are removed first.
	ixcolsrev = sort(ixcols, 'descend');
	
	ixcolsrev = ixcolsrev(:)';
	if all( ixcolsrev == (n0:-1:n0-p+1) )
		R = R(1:n0-p, 1:n0-p);
		Q = Q(:, 1:n0-p);
		return;
	end
	
	n = n0;
	for k = ixcolsrev
		if k == n
			R(:,end) = [];
			n = n-1;
			continue;
		end
		
		R(:, k) = [];
		for j = k : min(m,n-1)
			G = givens(R(j,j), R(j+1,j));
			R(j,j) =  G.r;
			R(j+1,j) = 0;
			
			Gblk = G.Block();
			R(j:j+1, j+1:n-1) = Gblk*R(j:j+1, j+1:n-1);

			if nargout > 1
				Q(:, j:j+1) = Q(:, j:j+1)*transpose(Gblk);
			end
		end
		n = n-1;
	end
	
	R = R(1:end-p, :);
	Q = Q(:, 1:end-p);
end
%%
