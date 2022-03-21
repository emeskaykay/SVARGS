% Ls_incremental::try_enter
%
% NOTE: If pivot_tolerance > 0, it is possible for some coefficients to be zeroed out.
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
function [This, ichoice, coeff, jzeroed] = try_enter(This, Xfull, bmask_try, pcutoff, maxentrycount)
ignore(maxentrycount)
	
	enforce_modifier(nargout > 0);
	
	if (nargin < 3) || isempty(bmask_try)
		bmask_try = ~(This.bmask_curr | This.bmask_exclude);
		
	else
		bmask_try = (bmask_try & ~(This.bmask_curr | This.bmask_exclude));
	end
	
	%%% Xtry = Xfull(:, bmask_try);
	m = nnz(bmask_try);
	%%%

	jzeroed = [];
	if (m == 0) %isempty(Xtry)
		%This.res_curr = This.ydata;
		%This.mse_curr = mean(This.res_curr.*This.res_curr);
		ichoice = 0;
		coeff = [];
		
	else
		N = size(This.ydata, 1);
		
		%%% m = size(Xtry, 2);
		%%%
		
		if isempty(This.Xcurr)
			qh = transpose(This.Xfull_sqrnorms(bmask_try));
			% Jcurr = zeros(1,N);
			Lcurr = zeros(0,N);
			
			dot_errtol = DP_ProbableErrorBound__(max(abs(This.ydata(:))), 0, 0, This.res_curr, 0, N);
			
		else
			depcount = 0;
			if This.pivot_tolerance > 0
				[~, ~, xperm, depcount] = qr_factor(This.Xcurr, This.pivot_tolerance);

				if depcount > 0
					xdep = xperm(end-depcount+1:end);
					jzeroed = This.ixcurr(xdep);
					This.bmask_curr(jzeroed) = false;
					This.Xcurr(:, xdep) = [];
					This.ixcurr(xdep) = [];
					This.coeffs_curr(xdep) = [];
					This.res_curr = This.ydata - This.Xcurr*This.coeffs_curr;
					This.mse_curr = sum(This.res_curr.*This.res_curr)/numel(This.res_curr);
				end
			end
			
			[Rcurr, Qcurr] = qr_factor(This.Xcurr, 0); % Unpermuted
			QcurrTx = transpose(Qcurr);
			% Jcurr = Rcurr\QcurrTx;
			
			if (depcount > 0) || isempty(This.Lprev)
				Lcurr = QcurrTx*Xfull;					% effectively dot product involved (1a)
				
			else
				qlast = QcurrTx(end, :);
				Lcurr = [This.Lprev; qlast*Xfull];	% effectively dot product involved (1b)
				
% 				assert(max(vecc(abs(Lcurr - QcurrTx*Xfull))) < 1e-10);
			end
			LTry = Lcurr(:, bmask_try);
			qh = This.Xfull_sqrnorms(bmask_try) - sum(LTry.*LTry, 1);	% effectively dot product involved (2)
			qh = transpose(qh);
			
			dot_errtol = DP_ProbableErrorBound__(max(abs(This.ydata(:))), Qcurr, LTry, This.res_curr, Rcurr, N);
		end
		
		%%% utry = transpose(This.res_curr)*Xtry;
		utry = transpose(This.res_curr)*Xfull;		% effectively dot product involved (3)
		utry = utry(bmask_try);
		%%%
		
		utry = transpose(utry);
		ds = (utry.*utry)./qh;
		ds(qh == 0) = 0;
		assert(all(isfinite(ds)));
		
		mse = This.mse_curr - ds/N;
		
		maxerr = 10*dot_errtol;
		mse( (mse <= 0) & (mse > -maxerr) ) = eps;
		if any(mse < 0)
			warning([newline, ...
				'Ls_incremental::try_enter, mse (= ', num2str(min(mse)), ...
				') is negative beyond tolerance of ', num2str(-maxerr), ...
				', probably due to numerical error.']);
			
			mse(mse < 0) = eps;
		end
		
		cxdelta = ones(m,1);
		direction = 1;
		[pvalues, gcstrengths] = model_pval(ones(m,1)*This.mse_curr, mse, cxdelta, N, direction);

		[pvalmin, jm] = min(pvalues);
		if nnz(pvalmin == pvalues) > 1
			[~, jm] = max(gcstrengths);
			pvalmin = pvalues(jm);
		end

		if isnan(pvalmin)
			ichoice = 0;
			coeff = 0;
			
		elseif pvalmin >= pcutoff
			ichoice = 0;
			coeff = [];

		else
			% Compute new coefficient.
			%coeff = utry(jm)/qh(jm);
			
			ixtries = find(bmask_try);
			ichoice = ixtries(jm);
			
			%%% xchoice = Xtry(:, jm);
			xchoice = Xfull(:, ichoice);
			%%%
			
			% Update regressor set.
			This.Xcurr(:, end+1) = xchoice;
			
			% Update existing coefficient
			% Small errors build up when updating coeffs. This can cause mse to
			% become non-positive, at which point we crash. So better to recompute.
 			% This.coeffs_curr = [This.coeffs_curr - coeff*Jcurr*xchoice; coeff]; % Update.
			[Rcurr, Qcurr] = qr_factor(This.Xcurr, 0);
 			This.coeffs_curr = Rcurr\(Qcurr' * This.ydata); % Recompute.
			coeff = This.coeffs_curr(end);
			
			% Update residuals, gc, pvalues, serial.
			This.res_curr = This.ydata - This.Xcurr*This.coeffs_curr;
			This.mse_curr = mse(jm);
			This.gc = gcstrengths(jm);
			This.pvalue = pvalmin;
			This.serial = This.serial + 1;
			
			% Update mask
			This.bmask_curr(ichoice) = true;
			This.ixcurr(end+1) = ichoice;
			
			% Update L0
			This.Lprev = Lcurr;
		end
	end
	
end
%%

% Probabilistic error bound in sequential dot product. With probability of lying
% outside K.sigma, the error in a.b is bounded by:
%		err <= K*(N+1)*w*eps, where w is the maximum absolute value of a and b combined.
% There are effectively two nested dot products in the code above DP_ProbableErrorBound__()
% Since X (=Xfull or Xcurr) is given data, we can assume that it is "correct". Also for
% simplicity we assume that R (=Rcurr) is correct upto rounding. Then the the error in
% Q (=Qcurr) is bounded by:
%		kap(R)*eps
% where kap(R) is the condition number, ||R||*||inv(R)||, of R.
% The error in c (=coeffs_curr) is then bounded by:
%		kap(R)*||inv(R)||*||y||*eps
% (since effectively c is contained in inv(R)*Q'*y - even though it is computed in a different way)
% This in turn implies that he error in res (=res_curr) is then bounded by:
%		kap(R)*||inv(R)||*||y||*||X||*eps (since res = y - X*c).
%	=	kap(R)^2 * ||y|| * eps (since Xcurr = Qcurr*R).
% Note that parallel error bound with binary tree adder is better (sqrt(N+1) instead
% of (N+1). See Huynh and Mucke, Error Analysis and Precision Estimation etc ...
function e = DP_ProbableErrorBound__(xmax, Q, L, res, R, N)
	K = 5;
	w1 = K*(N+1)*max( xmax, max(abs(Q(:))) );
	w2 = K*(N+1)*max( xmax, max(abs(L(:))) );
	w3 = K*(N+1)*max( xmax, max(abs(res)) );
	kap = cond(R);
	e = (w1*w2 + w3*w3*xmax*kap)*kap*eps('double');
end
%%
