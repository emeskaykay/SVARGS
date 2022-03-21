% class Svargs_params
% Parameters for SVARGS model fit.
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
classdef Svargs_params
	properties (SetAccess = 'private')
%		Number of components in the model
		k;
	
%		[ default = 1 : floor(sqrt(N)), where N is the length of the time series. ]:
%		A vector containing allowable lags >= 0.
		LagsAllowed;
	end
	
	properties
%		Statistical value to use for model comparison.
%		values: 'loglr' (log likelihood ratio) or 'gc' (granger causality)
		pType = 'gc';
		
%		[default: 0.2*N, where N is the number of time points]
%		Maximum number of coefficients allowed in any equation. If this is exceeded, the
%		fit is stopped. This can be set to a smaller number when the system is too large
%		to fit in a practical amount of time. By default it is set 0.2*N. Note that in any
%		particular equation, the value given to MaxCoeffsPerEquation is also bounded by
%		number of testable locations in that equation (due to the presence of exogeneous
%		variables or excluded locations). In general MaxCoeffsPerEquation must be small
%		in comparison to N, otherwise many of our assumptions are not valid.
		MaxCoeffsPerEquation = []; % If empty, will be set by fitter to 0.2*N.
		
%		Default pExit/pEntry ratio.
		ExitEntryRatio = 2;
		
%		Minimum rank to tolerate in QR factorization before omitting regressor(s).
		QrRankTol = 0;
		
%		[default: false]
%		If DetermineOrder is true, the fit function attempts to determine the model order
%		before starting the fit. It chooses the order from the set of lags specified in
%		LagsAllowed field. If DetermineOrder is false, the fit function uses all the lags
%		specified in LagsAllowed. In the latter case the order will be determined by the
%		as per the non-zero coefficients at the highest lag.
		DetermineOrder = false;

%		[default: false]
%		If true, a constant regressor is added to each equation fit. The output then contains
%		a non-zero intercept. This can also be a logical vector of length k indicating which
%		equations should include a constant regressor in the fit.
		bIncludeIntercept = false;
		
%		[default: @sbc_x]
%		Model quality measure. 
		hCriterion = @sbc_x;
		
% 		The Number of consecutive steps for which IC values must be greater than the current
% 		minimum IC value before the iteration stops.
		IcTailLength = 3;
		
% 		If non-zero, the iteration only stops if the last IcMonotoneTail consecutive IC values
%		are montone non-decreasing.
		IcMonotoneTail = 0;
		
%		If true, stores the history of coefficients tried.
		bHistory = false;
		
%		The maximum number of "Bad Loops" (Searches that become circular) before giving up.
%		If a "Bad Loop" is encountered, the code breaks out by discarding some (heuristically
%		chosen) coefficients. Bad Loops should not happen if pExit > pEntry.
		MaxBadLoops = 5;
	end
	
	properties (Dependent = true, SetAccess = 'private')
%		This is always first index in the TimeMask field. It is the Index of the
%		time series data that is to be treated as time zero. Data points prior to the
%		start index serve as initial data. If the start index is less than or equal to
%		the fit order (p), then the data is padded with (1 + p - StartIndex) leading
%		zeros. To prevent data padding, let StartIndex be at least p+1. We have observed
%		that data padding can have a very unpleasant effect on the accuracy of the
%		coefficients (particularly so when the coefficients are dense), so it is best
%		to leave this at the default value of p+1. However if the number of data points,
%		N, is small (say N < 10*p), then it may be better to set StartIndex=1. StartIndex
%		can be modified using SetStartIndex(). Modifying StartIndex will also modify
%		the TimeMask field.
		StartIndex;
		
%		[ default: [false(1,p), true(1,timelen-p+1))] ]
%		Boolean vector that indicates the locations of the regressor rows (in the set
%		of lagged regressors) that are to be used in the fit. The length of TimeMask
%		must equal the total number of time points in the data. It can be accessed
%		using GetTimeMask(). Setting the TimeMask can modify the StartIndex field.
		TimeMask;
		
%		Read only property. This is the number of lags in LagsAllowed.
		LagCount;
	end
	
	properties (Dependent = true)
%		[ default: None, ie: [] ]:
%		A logical (k x k x s) matrix, where s <= (Number of lags), and k is the number of
%		variables specifying locations that are added in an initial fit. Note that locations
%		present in bExclude are always removed from bInitial.
		bInitial;
		
%		[ default: None, ie: [] ]:
%		A logical (k x k x s) matrix, where s <= (Number of lags), and k is the number of
%		variables, specifying locations that will not be removed once they are added. If a
%		location is indicated both in bInitial as well as in bKeep, then that location will
%		always be fit with a coefficient regardless of its significance. Note however that
%		it is possible for the value of the coefficient to be 0.  Note that locations present
%		in bExclude are always removed from bKeep.
		bKeep;

%		[ default: None, ie: [] ]
%		A logical (k x k x s) matrix, where s <= (Number of lags), and k is the number of
%		variables, specifying which locations are to be excluded from the model. Note that
%		diagonal locations at lag zero are always excluded.
		bExclude;
	end
	
	properties (Access = 'private')
		btimemask__;
		tmaskset__ = false;
		b_initial__;
		b_keep__;
		b_exclude__;
	end
	
	properties %(Access = ?Svargs_fitter)
%		Lower bound on false detection rate. Default based on size of system.
		FdrMin;
		
%		Vector of estimates of the number of non-zero coefficients in each equation.
%		Changes along the fit. To be used only by fitting function.
		EstimatedCoeffCounts = [];
	end
	

	methods %(Access = ?Svargs_fitter)
		function ntl = GetNumberOfTestLocations(This, ieq)
			b = ~This.bExclude(ieq, :);
			ntl = nnz(b);
		end
		
		function This = SetEstimatedCoeffCount(This, ieq, M)
			enforce_modifier(nargout);
			This.EstimatedCoeffCounts(ieq) = M;
		end
	
		function fdcurr = GetFdrNext(This, fdcurr, ieq)
			enforce_modifier(nargout);
			fdmax = This.GetFdrMax(ieq);
			if fdcurr < fdmax
				fdrseq = This.GetFdrSequence(fdcurr, ieq);
				fdrseq = fdrseq(fdrseq > fdcurr);
				fdcurr = min(fdmax, fdrseq(1));
			end
		end
		
		% Upper bound on false detection rate. Varies along the fit.
		function fmax = GetFdrMax(This, ieq)
			% This ensures that pExit <= 0.5;
			ntl = This.GetNumberOfTestLocations(ieq);
			s = This.EstimatedCoeffCounts(ieq)/ntl;
			r = This.ExitEntryRatio;
			fmax = 1/(s*(r+1));
		end
		
		% Maximum p-value for adding a coefficient into the model.
		% Chosen based on FdrSeq during course of fit.
		function p = GetPentry(This, fdcurr, ieq)
			ntl = This.GetNumberOfTestLocations(ieq);
			s = This.EstimatedCoeffCounts(ieq)/ntl;
			r = This.ExitEntryRatio;
			p = fdcurr*s/(1 + s*(r-1)*fdcurr);
		end
		
		% Minimum p-value for removing a coefficient from the model.
		% Set to (ExitEntryRatio * pEntry).
		function p = GetPexit(This, fdcurr, ieq)
			p = This.ExitEntryRatio*This.GetPentry(fdcurr, ieq);
		end
	end
	
	methods (Access = 'private')
		function bArray = CheckIndicatorArray(This, bArray)
			q = This.LagCount;

			assert( islogical(bArray) && (ndims(bArray) <= 3) );
				
			% Resize bArray along the third dim.
			assert( size(bArray,1) == This.k );
			assert( size(bArray, 2) == This.k );
			assert( size(bArray,3) <= q );
			bArray = cat(3, bArray, false(This.k, This.k, q - size(bArray,3)) );
		end
		
		function n = NumFdrSteps(This, fdstart, ieq)
			n = ceil(1 + 2*log2(This.GetFdrMax(ieq)/fdstart));
		end
		
		% Sequence of FDR value to be tested.
		function fseq = GetFdrSequence(This, fdstart, ieq)
			n = This.NumFdrSteps(fdstart, ieq);
			fdmax = This.GetFdrMax(ieq);
			r = fdmax/fdstart;
			alpha = r^(1/n);
			fseq = fdstart*alpha.^(1:n);
			fseq(end) = fdmax;
		end
		
	end
	
	methods
		function This = Svargs_params(k, lags)
			assert( isscalar(k) && (k > 0) );
			assert( isnumeric(lags) && isvector(lags) && issorted(lags) && all(lags >= 0) );
			lags = sort(lags);
			
			This.k = k;
			This.LagsAllowed = lags;
			This.btimemask__ = [false(1, lags(end)), true];
			
			q = This.LagCount;
			This.b_initial__ = false(k,k,q);
			This.b_keep__ = false(k,k,q);
			This.b_exclude__ = false(k,k,q);
			if lags(1) == 0
				This.b_exclude__(:,:,1) = logical(eye(k));
			end
			
			This.FdrMin = 1/This.k/This.k/numel(This.LagsAllowed);
			This.FdrMin = min(0.1, This.FdrMin);
		end
		
		function params = Sub(This, ixvec, ixlagvec)
			if ~exist('ixlagvec', 'var'), ixlagvec = 1:This.LagCount; end
			assert( isvector(ixvec) && all(ixvec > 0) && all(ixvec <= This.k) && is_unique(ixvec) );
			assert( isvector(ixlagvec) && all(ixlagvec > 0) && all(ixlagvec <= This.LagCount) && is_unique(ixlagvec) );
			
			ixlagvec = sort(ixlagvec);
			
			params = This;
			params.k = length(ixvec);
			params.LagsAllowed = This.LagsAllowed(ixlagvec);
			
			params.bInitial = This.bInitial(ixvec, ixvec, ixlagvec);
			params.bKeep = This.bKeep(ixvec, ixvec, ixlagvec);
			params.bExclude = This.bExclude(ixvec, ixvec, ixlagvec);
		end
		
		% Accessors
		function b = GetTimeMask(This, ntimepts)
			if This.tmaskset__
				assert(ntimepts == numel(This.btimemask__), ...
					'ntimepts argument does not match length of previously set time mask.');
				b = This.btimemask__;
			else
				istart = This.StartIndex;
				b = [false(1, istart-1), true(1, ntimepts-istart+1)];
			end
		end
		
		function B = get.bInitial(This)
			B = This.b_initial__;
		end
		
		function B = get.bKeep(This)
			B = This.b_keep__;
		end
		
		function B = get.bExclude(This)
			B = This.b_exclude__;
		end
		
		function r = get.LagCount(This)
			r = numel(This.LagsAllowed);
		end
		
		function ix = get.StartIndex(This)
			ix = find(This.btimemask__, 1);
		end
		% ___________________
		
		function This = SetLagsAllowed(This, lags)
			enforce_modifier(nargout);
			
			assert( is_unique(lags) );
			assert( issorted(lags) );
			
			if ~all(ismember(lags, This.LagsAllowed))
				assert(isempty(This.b_initial__) || all(This.b_initial__(:) == 0));
				assert(isempty(This.b_keep__) || all(This.b_keep__(:) == 0));
				assert(isempty(This.b_exclude__) || all(This.b_exclude__(:) == 0));
				
				q = numel(lags);
				k_ = This.k;
				if ~isempty(This.b_initial__), This.b_initial__ = false(k_,k_,q); end
				if ~isempty(This.b_keep__), This.b_keep__ = false(k_,k_,q); end
				if ~isempty(This.b_exclude__), This.b_exclude__ = false(k_,k_,q); end
				
			else
				blagmask = ismember(This.LagsAllowed, lags);
				This.b_initial__ = This.b_initial__(:,:,blagmask);
				This.b_keep__ = This.b_keep__(:,:,blagmask);
				This.b_exclude__ = This.b_exclude__(:,:,blagmask);
			end
			This.LagsAllowed = lags;
			
			% Adjust the location of the first non-zero bit in TimeMask.
			This.btimemask__(1:lags(end)) = false;
			if length(This.btimemask__) == lags(end)
				This.btimemask__ = [This.btimemask__, true];
			end
		end
		
		function This = SetLagMultiplier(This, n)
			enforce_modifier(nargout);
			assert( (n >= 1) && (rem(n,1) == 0) );
			This.LagsAllowed = n*This.LagsAllowed;
			This.btimemask__ = vecr(This.btimemask__(ones(n,1), :));
		end
		
		function This = SetStartIndex(This, istart, bmask)
			enforce_modifier(nargout);
			assert(isscalar(istart) && isnumeric(istart));
			if nargin < 3
				bmask = true(size(This.btimemask__));
			end
			assert(all(size(bmask) == size(This.btimemask__)));
			This.btimemask__ = bmask;
			This.btimemask__(istart) = true;
			This.btimemask__(1:istart-1) = false;
		end
		
		function This = SetTimeMask(This, bmask)
			enforce_modifier(nargout);
			assert(isvector(bmask) && islogical(bmask));
			This.btimemask__ = vecr(bmask);
			This.tmaskset__ = true;
		end
		
		function This = set.bInitial(This, bx)
			This.b_initial__ = CheckIndicatorArray(This, bx);
			if ~isempty(This.b_exclude__)
				This.b_initial__ = (This.b_initial__ & ~This.b_exclude__);
			end
		end
		
		function This = set.bKeep(This, bk)
			This.b_keep__ = CheckIndicatorArray(This, bk);
			if ~isempty(This.b_exclude__)
				This.b_keep__ = (This.b_keep__ & ~This.b_exclude__);
			end
		end
		
		function This = set.bExclude(This, bx)
			This.b_exclude__ = CheckIndicatorArray(This, bx);
			if This.LagsAllowed(1) == 0
				This.b_exclude__(:,:,1) = (This.b_exclude__(:,:,1) | logical(eye(This.k)));
			end
			This.b_initial__ = (This.b_initial__ & ~This.b_exclude__);
			This.b_keep__ = (This.b_keep__ & ~This.b_exclude__);
		end
		
		function This = set.IcTailLength(This, n)
			assert(isnumeric(n) && isscalar(n) && (n > 0));
		end
		
		function This = set.IcMonotoneTail(This, n)
			assert(isnumeric(n) && isscalar(n) && (n >= 0));
		end
		
		function This = set.MaxBadLoops(This, m)
			assert( (m >= 1) && (rem(m,1) == 0) );
			This.MaxBadLoops = m;
		end
	end
	
end
%%
