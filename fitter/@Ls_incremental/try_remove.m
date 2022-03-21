% Ls_incremental::try_remove
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
function [This, ichoice] = try_remove(This, pcutoff, maxremcount)
ignore(maxremcount)

	enforce_modifier(nargout > 0);
	
	N = size(This.ydata, 1);
	
	[Rcurr, Qcurr] = qr_factor(This.Xcurr, 0);
	m = numel(This.ixcurr);
	msetry = zeros(m,1);
	% coeffs_curr_adjust = zeros(m-1, m);
	coeffs_curr_try = zeros(m-1, m);
	
	% Use QR updating for column deletion.
	for j = 1 : m
		% Do not remove regressors marked in bkeep.
		if This.bmask_keep(This.ixcurr(j))
			msetry(j) = Inf;
			continue;
		end
		
		[Rtry, Qtry] = qr_delete_cols(j, Rcurr, Qcurr);
		% [Qtry, Rtry] = qrdelete(Qcurr, Rcurr, j);
		QtryTx = transpose(Qtry);
		%Jtry = Rtry\QtryTx;

		x = This.Xcurr(:, j);
		c = This.coeffs_curr(j);
		% coeffs_curr_adjust(:, j) = c*Jtry*x;
		coeffs_curr_try(:,j) = Rtry\(QtryTx*This.ydata);
		
		Ltry = QtryTx*x;
		qh = sum(x.*x) - sum(Ltry.*Ltry);
		ds = c*c*qh/N;
		msetry(j) = This.mse_curr + ds;
	end
	
	cxdelta = ones(m,1);
	direction = -1;
	[pvals, gcstrengths] = model_pval(ones(m,1)*This.mse_curr, msetry, cxdelta, N, direction);
	
	[pvalmax, jm] = max(pvals);
	if nnz(pvalmax == pvals) > 1
		[~, jm] = min(gcstrengths);
		pvalmax = pvals(jm);
	end
	
	% Note: (jm == m) is not possible in theory, but can happen due to numerical error.
	% So test for it and do not remove if it occurs. In other words, we do not want
	% to remove the coefficient added at the previous step.
	if isempty(pvalmax) || (pvalmax <= pcutoff) || (jm == m)
		ichoice = 0;
		
	else
		% Update the object.
		% This.coeffs_curr = [This.coeffs_curr(1:jm-1, 1); This.coeffs_curr(jm+1:end, 1)] + coeffs_curr_adjust(:, jm);
		This.coeffs_curr = coeffs_curr_try(:,jm);
		This.Xcurr = [This.Xcurr(:, 1:jm-1), This.Xcurr(:, jm+1:end)];
		This.res_curr = This.ydata - This.Xcurr*(This.coeffs_curr);
		
		This.mse_curr = msetry(jm);
		This.gc = gcstrengths(jm);
		This.pvalue = pvalmax;
		This.serial = This.serial + 1;
		
		ichoice = This.ixcurr(jm);
		This.ixcurr(jm) = [];
		This.bmask_curr(ichoice) = false;
		
		% reset L0
		This.Lprev = zeros(0, N);
	end

end
%%

% Original code without using QR updating.
% 	for j = 1 : m
% 		Xtry = [This.Xcurr(:, 1:j-1), This.Xcurr(:, j+1:end)];
% 		
% 		This.pivot_tolerance = 0;
% 		[Rtry, Qtry, xperm, depcount] = qr_factor(Xtry, This.pivot_tolerance);
% 		
% 		Rtry = Rtry(1:end-depcount, 1:end-depcount);
% 		Qtry = Qtry(:, 1:end-depcount);
% 		xperm = xperm(1:end-depcount);
% 		QtryTx = transpose(Qtry);
% 		Jtry = Rtry\QtryTx;
% 		Jtry(xperm, :) = Jtry;
% 
% 		x = This.Xcurr(:, j);
% 		c = This.coeffs_curr(j);
% 		coeffs_curr_adjust(:, j) = c*Jtry*x;
% 		Ltry = QtryTx*x;
% 		qh = sum(x.*x) - sum(Ltry.*Ltry);
% 		ds = c*c*qh/N;
% 		msetry(j) = This.mse_curr + ds;
% 	end
