% SVARGS model fit.
%
% data:
%	(k x N) numeric matrix where k is the number of variables and N is the number of time points.
%
% fitparams:
%	An object of type Svargs_params. See Svargs_params.m
%
% dispobj [ default: 0 ]:
%	One of the following:
%	0:
%		No display output.
%
%	1 (true):
%		Stationary progress display on the screen.
%
%	[fid, b]:
%		A fileid obtained from a call to fopen, or 1 for std output.
%		All output is written to this file. If b is true, output is
%		detailed, otherwise each step is a single line. Default for
%		b is false.
%
%	gr:
%		A Graph object. Coefficient image is displayed on the graph during the fit.
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
classdef Svargs_fitter < handle
	properties (Access = 'private')
		% Constants
		reg_indicator_mask; % (k x lagcount) logical matrix indicating the regressors present in RegColumns.
		reg_indices; % Linear indices of regressor locations in regressor mask.
		testloc_count; % Number of coefficients to be tested.
		max_coeffcount;
		
		% Variable
		nbadloops = 0;
		fdr_curr;
		pentry_curr;
		dispinfo;
		eqn_orders;
	end
	
	properties (Access = 'private', Transient)
		hregcache;
	end
	
	properties (SetAccess = 'private')
		% Constants
		data;

		% Variable
		fitparams; % Some fields change temporarily during the course of the fit.
		stats;
		dispobj;
	end
	
	properties (Dependent = true, SetAccess = 'private')
		% Constant
		varcount;
	end
	
	methods
		function m = get.varcount(This)
			m = size(This.data, 1);
		end
	end
	
	methods (Static, Access = 'private')
		function b = has_turned__(ModHist, TailLength, MonotoneTail)
			[~, q0] = ModHist.GetBestLevel();
			[~, qvec] = ModHist.IcTail();
			b = (numel(qvec) >= TailLength) && all(qvec >= q0);
			m = min(MonotoneTail, TailLength);
			b = b && all(diff(qvec(end-m:end)) >= 0);
		end
	end
	
	methods (Access = 'private')
		% Currently implemented for block size = 1 only.
		function [Y, RegColumns, bExcludeIntercept] = Initialize(This)
			btmask = This.fitparams.GetTimeMask(size(This.data, 2));
			
			This.reg_indicator_mask = any(false == This.fitparams.bExclude(:,:), 1);
			This.reg_indicator_mask = reshape(This.reg_indicator_mask, This.fitparams.k, This.fitparams.LagCount);
			
			[RegColumns, Y] = stack_regressors(This.data, This.hregcache, btmask, This.fitparams.LagsAllowed, This.reg_indicator_mask, false);
			This.reg_indices = find(This.reg_indicator_mask);
			This.testloc_count = numel(This.reg_indices);
			
			N = size(Y, 1);
			RegColumns = [ones(N,1), RegColumns];
			bExcludeIntercept = vecc(false == This.fitparams.bIncludeIntercept);
			if isscalar(bExcludeIntercept)
				bExcludeIntercept = bExcludeIntercept(ones(This.fitparams.k, 1), 1);
			end
			
			Nused = size(Y,1); % Nused may be different from N if there is a non-trivial time mask or istart > 1.
			This.stats = Svargs_hist(Nused, This.fitparams.hCriterion);
			
			% In theory the maximum number of coefficients, m, in any equation must
			% satisfy m << N, otherwise many of our assumptions are not valid.
			if isempty(This.fitparams.MaxCoeffsPerEquation)
				This.fitparams.MaxCoeffsPerEquation = (N-1)/5;
			end
			for ieq = 1 : This.fitparams.k
				This.max_coeffcount(ieq) = min(This.fitparams.MaxCoeffsPerEquation, This.fitparams.GetNumberOfTestLocations(ieq));
				This.fitparams = This.fitparams.SetEstimatedCoeffCount(ieq, 1); % Initialize
			end
			
			This.eqn_orders = zeros(This.fitparams.k, 1);
			This.dispinfo = initialize_display(This.dispobj, This.fitparams);
		end
		
		function [Coeffs, Res, Sigma, Intercept] = Finalize(This)
			k = This.fitparams.k;
			Lags = This.fitparams.LagsAllowed;
			q = numel(Lags);
			
			Cmat = This.stats.GetBestCoeffs();
			Intercept = Cmat(:,1);
			Cmat = Cmat(:,2:end);
			Coeffs = zeros(k,k,q);
			Coeffs(:, This.reg_indicator_mask) = Cmat;
			
			[RegColumns1, XY] = stack_regressors(This.data, This.hregcache, 1, Lags, 'all', false);
			Res = XY - RegColumns1*transpose(Coeffs(:,:));
			Res = transpose(Res);
			Sigma = covariance(Res);
			
			This.dispinfo.order = max(This.eqn_orders);
			finalize_display(This.dispinfo, [], Coeffs);
		end
		
		% Test the coefficients that are are not present in the current model.
		function [LsFit, bEntered, bZeroed, coeff] = TryEnter__(This, RegColumns, bMaskTry, LsFit, ModHist)
			ieq = ModHist.objid;
			[LsFit, jchoice, coeff, jzeroed] = LsFit.try_enter(RegColumns, bMaskTry, This.fitparams.GetPentry(This.fdr_curr, ieq), 1);

			bZeroed = false;
			if (jchoice == 0)
				bEntered = false;

			elseif coeff == 0
				disp('Entered coefficient == 0. Likely cause: One of the response variables has all zeros. ');
				bEntered = false;

			else
				enter__(This, jchoice, coeff, LsFit, ModHist)
				bEntered = true;
				
				for j0 = jzeroed
					remove__(This, j0, LsFit, ModHist);
					bZeroed = true;
				end
			end
		end
		
		function enter__(This, jchoice, coeff, LsFit, ModHist)
			k = This.fitparams.k;
			
			if jchoice == 1 % first column is constant regressor.
				kfrom = 0;
				lag = 0;
				
			else
				ixreg = This.reg_indices(jchoice - 1);

				kfrom = rem( ixreg-1, k ) + 1;
				all_lags = This.fitparams.LagsAllowed;
				lag = all_lags(floor((ixreg-1)/k) + 1);
			end

			id.to = ModHist.objid;
			id.from = kfrom;
			id.lag = lag;

			% Add new coefficient and update existing ones.
			ModHist.EnterCoeffs(id, jchoice, coeff, vecr(LsFit.ixcurr), vecr(LsFit.coeffs_curr), LsFit.mse_curr);
		end
		
		% Test the coefficients that are are present in the current model for removal
		function [LsFit, bRemoved] = TryRemove__(This, LsFit, ModHist)
			ieq = ModHist.objid;
			[LsFit, jchoice] = LsFit.try_remove(This.fitparams.GetPexit(This.fdr_curr, ieq), 1);

			if jchoice == 0
				bRemoved = false;

			else
				remove__(This, jchoice, LsFit, ModHist)
				bRemoved = true;
			end
		end
		
		function remove__(This, jchoice, LsFit, ModHist)
			k = This.fitparams.k;
			
			if jchoice == 1 % first column is constant regressor.
				kfrom = 0;
				lag = 0;
				
			else
				ixreg = This.reg_indices(jchoice - 1);

				kfrom = rem( ixreg-1, k ) + 1;
				all_lags = This.fitparams.LagsAllowed;
				lag = all_lags(floor((ixreg-1)/k) + 1);
			end

			id.to = ModHist.objid;
			id.from = kfrom;
			id.lag = lag;
			
			% Remove coefficient and update existing ones.
			ModHist.EnterCoeffs(id, jchoice, 0, vecr(LsFit.ixcurr), vecr(LsFit.coeffs_curr), LsFit.mse_curr);
		end
		
		% Stores history of differences and checks for circularity.
		function [ClocsHist, Circ] = CheckCircular__(This, Clocs, ClocsHist)
			nb = numel(Clocs);
			if isempty(ClocsHist)
				ClocsHist = sparse(nb, 0);
			end
			ClocsOld = sparse(nb, 1); % zeros.
			ClocsNew = sparse( double(Clocs(:)) );

			i = 1;
			Circ = false;
			nSteps = size(ClocsHist,2);
			while (Circ == false) && (i <= nSteps)
				ClocsOld = ClocsOld + ClocsHist(:,i); % Cblocks value at step i.
				Circ = all(ClocsNew == ClocsOld);
				i = i + 1;
			end

			if Circ == false
				% Store differences only.
				ClocsHist = [ClocsHist, ClocsNew - ClocsOld];

			elseif This.nbadloops < This.fitparams.MaxBadLoops
				ClocsHist = [];
				This.nbadloops = This.nbadloops + 1;
				This.fitparams.EntryExitRatio = This.fitparams.EntryExitRatio + 1;

			else
				print_obj(1, 'FitVARBlockwise: Search has become circular. ');
				print_obj(1, 'Try Lower ''pEntry'' and higher ''pExit'' values. \n');
				error('FitVARBlockwise:Circular search. ');
			end
		end
		
		function SetDisplayInfo(This, LsFit, ModHist, ieq)
			if isfield(This.dispinfo, 'G')
				This.dispinfo.ImgInfo = This.dispinfo.ImgInfo.setData(ModHist.coeffs);
				This.dispinfo.G.PlotImageUpdate(This.dispinfo.ImgInfo, 1);
			end

			This.dispinfo.serial = ModHist.serial;

			This.dispinfo.eq = 0;
			This.dispinfo.lag = NaN;
			if ~isempty(ModHist.entryid)
				This.dispinfo.eq = ModHist.entryid.to;
				This.dispinfo.lag = ModHist.entryid.lag;
				This.dispinfo.from = ModHist.entryid.from;
			end

			This.dispinfo.coeff_count = ModHist.coeff_count;
			This.dispinfo.order = This.eqn_orders(ieq);
			This.dispinfo.ic_taillength = ModHist.IcTailLength();
			This.dispinfo.fdrcurr = This.fdr_curr;
			This.dispinfo.pentrycurr = This.pentry_curr;

			if isfield(This.dispinfo, 'ShowDetail') && This.dispinfo.ShowDetail
				This.dispinfo.pvalue = LsFit.pvalue;
				This.dispinfo.mse = LsFit.mse_curr;
				This.dispinfo.gc = LsFit.gc;

				if ModHist.cxdelta < 0
					This.dispinfo.action = 'REMOVED';

				elseif ModHist.cxdelta > 0
					This.dispinfo.action = 'ADDED';

				else
					assert(false);
				end
			end
		end
		
	end
	
	methods
		function This = Svargs_fitter(dx, fitprms, disp)
			% Set default arguments.
			if (nargin < 3), disp = 1; end
			
			assert( isa(fitprms, 'Svargs_params') );
			
			This.hregcache = VariablesCache('Regressors', 'RegressorLags');
			This.data = dx;
			This.fitparams = fitprms;
			This.dispobj = disp;
		end

		function [Coeffs, Res, Sigma, Intercept] = Run(This)
			[Y, RegColumns, bExcludeIntercept] = This.Initialize();
			qmax = numel(This.fitparams.LagsAllowed);

			R = qr_factor(RegColumns);
			RegColumnSqrNorms = sum(R.*R, 1);

			R = qr_factor(Y);
			Ysqrnorms = sum(R.*R, 1);
			
			FitParamsOrig = This.fitparams;
			for ieq = 1 : This.varcount
				binitial = [false, This.fitparams.bInitial(ieq, This.reg_indicator_mask)];
				bkeep = [false, This.fitparams.bKeep(ieq, This.reg_indicator_mask)];
				bexclude = [bExcludeIntercept(ieq), This.fitparams.bExclude(ieq, This.reg_indicator_mask)];
				
				LsFit = Ls_incremental(Y(:,ieq), Ysqrnorms(ieq), RegColumnSqrNorms, bkeep, bexclude);
				LsFit.pivot_tolerance = This.fitparams.QrRankTol;
				LsFit = LsFit.initialize(RegColumns, binitial);
				
				coeffs_initial = zeros(size(LsFit.bmask_curr));
				coeffs_initial(LsFit.bmask_curr) = LsFit.coeffs_curr;
				ModHist = This.stats.NewEquationHist(coeffs_initial, LsFit.mse_curr);
				
				if all(This.reg_indicator_mask(ieq, :) == false)
					continue;
				end
				
				ClocsHist = [];
				
				if This.fitparams.DetermineOrder
					NzLagCount = 1; % Start with single lag and gradually increase.
				else
					NzLagCount = qmax; % Standard option
				end
				
				% Enter/remove coefficients at each order until model remains unchanged
				% Then increment order. Repeat until all specified orders are covered.
				bEntered = true;
				bRemoved = true;
				bChanged = (bEntered || bRemoved);
				bTurned = false;
				This.fdr_curr = This.fitparams.FdrMin;
				This.pentry_curr = This.fitparams.GetPentry(This.fdr_curr, ieq);
				while bChanged && ~bTurned
					% Try to add and display statistics if entered.
					bMaskTry = (LsFit.bmask_curr == false);
					[LsFit, bEntered, bZeroed, coeff_entry] = This.TryEnter__(RegColumns, bMaskTry, LsFit, ModHist);
					if (bEntered || bZeroed)
						This.eqn_orders(ieq) = max(This.eqn_orders(ieq), ModHist.entryid.lag);
						This.SetDisplayInfo(LsFit, ModHist, ieq);
						This.dispinfo.last_entry = coeff_entry;
						update_display(This.dispinfo);
					end
					
					% Try to remove and display statistics if removed.
					[LsFit, bRemoved] = This.TryRemove__(LsFit, ModHist);
					if (bRemoved)
						This.dispinfo.last_entry = 0;
						This.eqn_orders(ieq) = max(This.eqn_orders(ieq), ModHist.entryid.lag);
						This.SetDisplayInfo(LsFit, ModHist, ieq);
						update_display(This.dispinfo);
					end

					if This.fitparams.bHistory
						ClocsHist = This.CheckCircular__(ModHist.coeffs ~= 0, ClocsHist);
					end

					bChanged = (bEntered || bZeroed || bRemoved);
					bTurned = Svargs_fitter.has_turned__(ModHist, This.fitparams.IcTailLength, This.fitparams.IcMonotoneTail);
					
					if ~bChanged
						% If no coefficients were added or removed:
						%	1. Refine sparsity estimate (which will affect the pEntry and pExit values).
						%	2. Increase the lag count. If lag count is already maximum, increase the FDR value.
						This.fitparams = This.fitparams.SetEstimatedCoeffCount(ieq, max(1, LsFit.curr_count()));

						if NzLagCount < qmax
							NzLagCount = NzLagCount + 1;
							bChanged = true; % Order has changed.

 						elseif Svargs_fitter.has_turned__(ModHist, This.fitparams.IcTailLength, This.fitparams.IcMonotoneTail)
 							bTurned = true;
							
						elseif (This.fdr_curr < This.fitparams.GetFdrMax(ieq)) ...
								&& (LsFit.curr_count() < This.max_coeffcount(ieq))
							This.fdr_curr = This.fitparams.GetFdrNext(This.fdr_curr, ieq);
							This.pentry_curr = This.fitparams.GetPentry(This.fdr_curr, ieq);
							bChanged = true; % Need this here because pEntry and pExit have now changed.
						end
					end
				end
				bforce = true;
				This.SetDisplayInfo(LsFit, ModHist, ieq);
				update_display(This.dispinfo, bforce);
				
				This.fitparams = FitParamsOrig;
				neweq_display(This.dispinfo);
			end

			[Coeffs, Res, Sigma, Intercept] = Finalize(This);
		end
		
	end
	
end
%%
