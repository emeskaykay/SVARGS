% Ls_incremental
% Performs an incremental least squares regression fit to the data.
% Statistical information about the changes are returned.
%
% Usage:
%	Obj = Ls_incremental(y, ynorm, regressor_norms)
%
% y:
%	Vector data (Response vector)
%
% ynorm:
%	Euclidean norm of y.
%
% regressor_norms:
%	Row vector of squared norms of the full set of regressors.
% _______________________________________________________________________
% Once the object is initialized, incremental fits can be performed using one of:
%			[Obj IsEntered] = Obj.try_enter(Xfull, bmask_try, pcutoff, maxentrycount);
% OR		[Obj IsRemoved] = Obj.try_remove(pcutoff, maxremcount);
%
% A call to Obj.try_enter(Xfull, bmask_try, pcutoff, maxentrycount) will try each
% regressor in Xfull(:, bmask_try) and enter the maxentrycount most significant
% coefficients whose pvalues are less than pcutoff. By default maxentrycount = 1.
% 
% A call to try_remove(pcutoff, maxremcount) will remove the maxrecount least
% significant coefficients whose pvalues are greater than pcutoff. By default
% maxremcount = 1.
% 
% try_enter() or try_remove() can be called any number of times. At each step the
% statistics can be examined. See try_enter.m and try_remove.m.
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
classdef Ls_incremental
	
	properties (Access = 'private')
		Xfull_sqrnorms;
		Lprev;
		Xcurr;
	end
	
	properties (SetAccess = 'private')
		ydata;			% vector data
		
		bmask_keep;		% Mask indicating non-removable regressors with respect to full regressor set.
		bmask_curr;		% Mask indicating current regressors (including non-removable ones) with respect to full regressor set (Xfull)
		bmask_exclude; % Mask indicating non-enterable regressors with respect to full regressor set.
		ixcurr;			% ixcurr are the indices, in order, of the non-zero elements of bmask_curr.
		coeffs_curr;	% Vector of coefficients corresponding to each of the regressors in Xfull(bmask_curr).
		res_curr;		% Vector of residuals.
		mse_curr;		% Mean square of residuals.
		
		gc;				% GC corresponding to the entered/removed coefficient.
		pvalue;			% pvalue for the entered/removed coefficient.
		
		serial;
	end
	
	properties
		pivot_tolerance = 0;
	end
	
	methods (Access = 'public')
		% Constructor.
		function This = Ls_incremental(y, ysqrnorm, regnorms, bkeep, bexclude)
			assert(isnumeric(y) && iscolumn(y));
			N = numel(y);
			
			assert(isnumeric(regnorms) && isrow(regnorms));
			K = numel(regnorms);
			
			if (nargin < 5) || isempty(bexclude), bexclude = false(1, K); end
			if (nargin < 4) || isempty(bkeep), bkeep = false(1, K); end
			
			assert(isrow(bkeep) && islogical(bkeep));
			assert(isrow(bexclude) && islogical(bexclude));
			assert(numel(bkeep) == K);
			assert(numel(bexclude) == K);
			assert(~any(bkeep & bexclude));
			assert(isrow(regnorms));
			
			This.Xfull_sqrnorms = regnorms;
			This.Lprev = zeros(0, N);
			This.ydata = y;
			
			This.Xcurr = zeros(N, 0);

			This.bmask_curr = false(1, K);
			This.bmask_keep = bkeep;
			This.bmask_exclude = (bexclude | (This.Xfull_sqrnorms == 0));
			This.ixcurr = [];

			This.coeffs_curr = zeros(0,1);
			This.res_curr = y;
			This.mse_curr = ysqrnorm/N;
			
			This.serial = 0;
		end

		function This = initialize(This, Xfull, binitial)
			enforce_modifier(nargout > 0);
			
			K = numel(This.Xfull_sqrnorms);
			assert(size(Xfull,2) == K);
			
			assert(isrow(binitial) && islogical(binitial));
			assert(numel(binitial) == K);
			assert(~any(binitial & This.bmask_exclude));
			
			if any(binitial)
				This.Xcurr = Xfull(:, binitial);
				This.ixcurr = find(binitial);
				
				[This.coeffs_curr, ~, ~, ixdep] = qr_regress(This.ydata, This.Xcurr, [], false, This.pivot_tolerance);
				
				This.ixcurr(ixdep) = [];
				This.bmask_curr(This.ixcurr) = true;
				This.Xcurr(:, ixdep) = [];
				This.coeffs_curr(ixdep) = [];
				
				This.coeffs_curr = vecc(This.coeffs_curr);
				This.res_curr = (This.ydata - This.Xcurr*This.coeffs_curr);
				This.mse_curr = sum(This.res_curr.*This.res_curr)/numel(This.res_curr);
			end
		end
		
		% Incremental fit.
		[This, ichoice, coeff, jzeroed] = try_enter(This, Xfull, bmask_try, pcutoff, maxentrycount)
		[This, ichoice] = try_remove(This, pcutoff, maxremcount)
		
		function m = curr_count(This)
			m = numel(This.ixcurr);
		end
		
		function This = clear_cache(This)
			This.Lprev = [];
			This.Xfull_sqrnorms = [];
			This.Xcurr = [];
		end
	end
	
end
%%
