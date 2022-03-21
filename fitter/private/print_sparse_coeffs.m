% print_sparse_coeffs(File, C, lags)
%
%
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
function print_sparse_coeffs(File, C, lags)
	k = size(C,1);
	q = numel(lags);
	
	for i = 1 : q
		if ~any(C(:,:,i)), continue; end
		
		print_obj(File, ['At lag ', num2str(lags(i)), ': ']);
		
		for r = 1 : k
			if ~any(C(:,r,i)), continue; end
			
			for s = 1 : k
				if C(s,r,i) == 0, continue; end
				
				coeffstr = ['b', num2str(s), '_', num2str(r)];
				str = [ '\t ', num2str(r, '% 4u'), ' ---> ', num2str(s, '% 4u'), ...
					' (', coeffstr, ' = ', num2str(C(s,r,i), '% 6.2f'), ')' ];
				print_obj(File, str);
			end
		end
		print_obj(File);
	end
	
	print_obj(File, ['\t Total coefficient count: \t', num2str(nnz(C))]);
	
end
%%
