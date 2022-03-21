% dispinfo = initialize_display(dispflag, FitParams)
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
function dispinfo = initialize_display(dispflag, fitparams)

	dispinfo.ic_criterion = func2str(fitparams.hCriterion);
	dispinfo.y_lags = fitparams.LagsAllowed;
	
	if ~isempty(dispflag) && isnumeric(dispflag) && ~all(dispflag == 0)
		if isscalar(dispflag) && dispflag
			% Stationary Progress_display for screen only.
			dispinfo.Progress = Progress_display('SVARGS fit', ...
				{'step', 'equation', 'p-entry', 'order', 'coeff-count', 'ictail'}, 'min_update_interval', 0.5);
			
		else
			File = dispflag(1);
			print_obj(File, '____________Begin SVARGS Fit Details____________');

			lags = fitparams.LagsAllowed;
			
			c = ['0', '1'];
			coeffcount = nnz(fitparams.bInitial);

			print_obj(File, '_________________Initial Model ___________________ ');
			print_obj(File, ['CoefficientCount: ', num2str(coeffcount)] );
			print_obj(File, ['Lags: ', num2str(lags)] );
			print_obj(File, ['Exogeneous Variable Count: ', num2str(nnz(fitparams.bExogeneous))] );
			print_obj(File, ['Exogeneous Mask: ', c(1+fitparams.bExogeneous)] );
			print_obj(File);

			StatValueType = 'Granger Causality';

			print_obj(File, '________________p-Value Cutoff''s_________________ ');
			print_obj(File, ['** p-value for entry \t = ', num2str(fitparams.pEntry)]);
			print_obj(File, ['** p-value for removal \t = ', num2str(fitparams.pExit)]);
			print_obj(File, '__________________________________________________');
			print_obj(File, ['p-value will be based on ', StatValueType]);
			print_obj(File);

			dispinfo.F = File;
			dispinfo.ShowDetail = logical(dispflag(2));
		end
	end
	
end
%%
