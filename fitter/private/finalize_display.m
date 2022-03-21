% finalize_display(dispinfo, CX, CY)
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
function finalize_display(dispinfo, CX, CY)
	if ~isempty(dispinfo)
		if isfield(dispinfo, 'Progress')
			dispinfo.Progress.complete();
			
			msgstr = ['Result after using the ''', dispinfo.ic_criterion, ''' criterion: ', ...
				' Total coefficient count: ', num2str(nnz(CX) + nnz(CY)), ', Order: ', num2str(dispinfo.order)];
			disp(msgstr);
			
		elseif isfield(dispinfo, 'F')
			File = dispinfo.F;

			% Display each coefficient.
			print_obj(File, '_________________Coefficients_____________________');
			print_sparse_coeffs(File, CX, 0);
			print_sparse_coeffs(File, CY, dispinfo.y_lags);
			print_obj(File);

			print_obj(File, '||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||');
			print_obj(File);
		end
	end
	
end
%%
