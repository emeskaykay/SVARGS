% update_display(dispinfo, bforce)
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
function update_display(dispinfo, bforce)

	if (nargin < 2) || isempty(bforce), bforce = false; end
	
	if ~isempty(dispinfo)
		if isfield(dispinfo, 'Progress')
			nvec = [dispinfo.serial, dispinfo.eq, dispinfo.pentrycurr, dispinfo.order, dispinfo.coeff_count, dispinfo.ic_taillength];
			dispinfo.Progress.update(nvec, bforce);
			
		elseif isfield(dispinfo, 'F')
			File = dispinfo.F;
			
			print_obj(File, ['Step # ', num2str(dispinfo.serial), ...
				', Eq. # ', num2str(dispinfo.eq), ...
				', Lag ', num2str(dispinfo.lag), ...
				', Current Fdr ', num2str(dispinfo.fdrcurr), ...
				', Coeff count: ', num2str(dispinfo.coeff_count)]);

			if dispinfo.ShowDetail
				print_obj(File, [action, ' interaction at lag ', num2str(dispinfo.lag)]);
				print_obj(File, ['\t From variable ', num2str(dispinfo.from)]);
				print_obj(File, ['\t   To variable ', num2str(dispinfo.eq)]);
				print_obj(File, ['\t p-value ', num2str(dispinfo.pvalue)]);
				if ~isempty(BlockSpec), print_obj(File, ['\t Coefficient = ', num2str(dispinfo.last_entry, 5)]); end
				print_obj(File);
				print_obj(File, ['\t Error estimate = ', num2str(dispinfo.mse)]);
				print_obj(File, ['\t ', 'Granger Causality = \t', num2str(dispinfo.gc)]);
				print_obj(File);
			end
		end
		
	end
	
end
%%
