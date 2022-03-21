% BlockCompanion
% Returns the corresponding AR(1) matrix of an AR(p) process.
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
function Mat = BlockCompanion(C)

	k = size(C,1);
	
	if isempty(C)
		Mat = zeros(k,k,0);
		return;
	end
	
	p = size(C, 3);
	Mat = [ C(:,:); [eye(k*(p-1)), zeros(k*(p-1), k)] ];
	
end
%%
