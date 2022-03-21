% BlockCompanionSparse
% Returns the corresponding AR(1) matrix of an AR(p) process as a sparse matrix.
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
function Mat = BlockCompanionSparse(Phi)
	k = size(Phi,1);
	
	if isempty(Phi)
		Mat = sparse(zeros(k,k,0));
		return;
	end

	p = size(Phi,3);
	Mat = sparse(Phi(:,:));
	I0 = speye(k*p - k);
	Z0 = sparse([], [], [], k*p - k, k, 0); % sparse ((k*p - k) x k) matrix of zeros.
	Mat = [ Mat; [I0, Z0] ];
end
%%
