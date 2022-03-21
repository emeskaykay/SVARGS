% [sr, checkfailed] = SpectralRadius(Op)
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
function [sr, checkfailed] = SpectralRadius(Op)
	
	p = Op.Order();
	
	C = Op.Coeffs(1:p, 'reduced', true);
	minsparsity = 50;
	eig_sp = (nnz(C) < numel(C)/minsparsity);
	
	if eig_sp
		A = BlockCompanionSparse(C);
	else
		A = BlockCompanion(C);
	end
	
	sr = Inf;
	checkfailed = false;
	if eig_sp
		opts.maxit = 1000;
		opts.disp = 0;
		try
			[~, l, flag] = eigs(A, 1, 'lm', opts);
			
			if any(isnan(l))
				print_x('%s', 'eigs() returned NaN value.');
				eig_sp = false;
			elseif flag ~= 0
				eig_sp = false;
			else
				sr = abs(l);
			end
	
		catch err
			print_x('%s\n', err.message);
			eig_sp = false;
		end
	end
	
	if (eig_sp == false)
		print_x('%s\n', 'Trying eig() ... ');
		try
			l = Evalues(Op);
			sr = max(abs(l));
			
		catch err
			print_x('\n');
			msg = sprintf('%s\n', err.message);
			msg = [msg, sprintf('%s\n', 'Failed to compute spectral radius. ')];
			disp_x(msg);
			checkfailed = true;
			return;
		end
	end
end
%%
