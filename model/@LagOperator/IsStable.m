% IsStable
% Usage: IsStable()
% Checks whether the system is stable by determining if powers of the lag
% operator converge or diverge.
%
% MyOne [default: 1 - eps('single')]:
%	
% bUseHeuristic [default: false]:
%	If true, heuristics can be used to determine instability. If bUseHeuristic == true, and
%	the returned value of 'b' is is zero (unstable), there is a small chance that the operator
%	is actually stable. If set to true, the computation should be much faster.
%
% bSuppressUnstable [default: false]:
%	Suppress display if result is unstable. This can be useful when carrying out a large
%	number of trials in order to find a stable model.
%
% Returns:
% b:
%	Boolean value indicating if the model is stable.
% If b == true:
%	The model is guaranteed to be stable.
% If b == false:
%	If bInfinite and bHeuristic are both false, the model is guaranteed to be unstable.
%	If any of bInfinite or bHeuristic is true, the model is unstable with high probability.
%
% bInfinite:
%	The second boolean return value bInfinite is false unless the power iteration reached
%	infinity (Inf) in which case bInfinite is true. If bInfinity == true, the returned value
%	of 'b' is is false (unstable) but there is an (extremely) small chance that the operator
%	is actually stable.
%
% bHeuristic:
%	True if heuristic was actually used.
%
% loopexp:
%	The last power of A that was checked.
%
% sr:
%	Spectral radius.
%	If the stability/instability is confirmed by an eigenvalue check then sr is the spectral
%	radius of the model. If the eigenvalue check is indecisive, then sr is set to 1. If the
%	model is unstable, then sr is only computed if the sr output argument is present. So it
%	is much faster to omit the sr in the output arguments if it is not required.
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
function [b, bInfinite, bHeuristic, loopexp, sr] = IsStable(Op, MyOne, bUseHeuristic, bSuppressUnstable)
NCRIT = 300;

	if (nargin < 4) || isempty(bSuppressUnstable), bSuppressUnstable = false; end
	if (nargin < 3) || isempty(bUseHeuristic), bUseHeuristic = false; end
	if (nargin < 2) || isempty(MyOne), MyOne = 1 - eps('single'); end
	
	
	if IsEmpty(Op)
		b = true;
		bInfinite = false;
		bHeuristic = false;
		loopexp = 0;
		sr = 0;
		return;
	end
	
	dataid = 'stable_indicator';
	data = Op.hCache.Get(dataid);
	if ~isempty(data) && ((nargout <= 4) || isempty(data.sr))
		b = data.b;
		sr = data.sr;
		bInfinite = data.bInfinite;
		bHeuristic = data.bHeuristic;
		loopexp = data.loopexp;
		return;
	end
	
	bInfinite = false;
	bHeuristic = false;
	loopexp = 0;
	sr = [];
	
	k = NumVars(Op);
	p = Order(Op);
	q = length(Op.lags);
	n = k*p;
	checkfailed = false;
	
	if (nargout > 4) % If spectral radius is desired, might as well try to compute it to determine stability.
		[b, checkfailed, sr] = CheckSpectralRadius__(Op);
		if checkfailed
			PrintResult__('LagOperator:IsStable: failed to compute spectral radius. ');
			PrintResult__('Attempting to determine stability: ', '%s');
			
		else
			if b
				PrintResult__(['Stable, Spectral Radius < 1 (', num2str(sr), ')']);
			else
				PrintResult__(['Unstable, Spectral Radius >= 1 (', num2str(sr), ')'], [], bSuppressUnstable);
			end
			return;
		end
	end

	minsparsity = 5*feature('numcores');
	C = Op.Coeffs(1:p, 'reduced', true);
	if (n > NCRIT) && (minsparsity*Op.CoeffCount() < numel(C))
		A = BlockCompanionSparse(C);
	else
		A = BlockCompanion(C);
	end
	
	if ~checkfailed && ~issparse(A)
		% Power iteration for non-sparse A can lead to large errors, so just use eig().
		l = eig(A);
		sr = max(abs(l));
		b = (sr < 1);

% 	elseif ~MemoryCheckMsg(2*n*k*q)
% 		error('LagOperator::IsStable: Execution cancelled. ');

	else
		[b, bInfinite, bHeuristic, loopexp] = PowerIteration__(A, k, minsparsity, MyOne, bUseHeuristic, bSuppressUnstable);
	end
	
	% If stable, do a confirmatory check
	if ~checkfailed && b
		[b, checkfail, sr] = CheckSpectralRadius__(Op);
		if checkfail
			PrintResult__('LagOperator:IsStable: Could not confirm, but trust norm. Assumed Stable. ');
			sr = 1;
			b = true;
			
		elseif b
			PrintResult__('OK');
			
		else
			PrintResult__('But check gives Spectral Radius >= 1). Assumed Unstable', [], bSuppressUnstable);
		end
	end
	
	data.b = b;
	data.sr = sr;
	data.bInfinite = bInfinite;
	data.bHeuristic = bHeuristic;
	data.loopexp = loopexp;
	Op.hCache.Set(dataid, data);
	
end

function [b, bInfinite, bHeuristic, loopexp] = PowerIteration__(A, k, minsparsity, MyOne, bUseHeuristic, bSuppressUnstable)

	b = false;
	bInfinite = false;
	bHeuristic = false;
	loopexp = 1;
	
	MyInf = Inf;
	n = length(A);
	k_upper = k;
	tr_A2 = trace(A);

	while (1)
		if (k_upper <= n/5)
			A = SquareVAR1TypeMatrix__(A, k_upper);
		else
			A = A*A;
		end
		k_upper = 2*k_upper;

		if (minsparsity*nnz(A) >= numel(A))
			% Change to full form.
			A = full(A);
		end

		a1 = norm(A, 1);
		ainf = norm(A, Inf);
		if min(a1, ainf) < MyOne
			PrintResult__('Stable (1 or inf norm converges). ');
			b = true;
			break;
		end

		% Very important since trace may not diverge.
		% (Inf - Inf) can lead to NaN, so check for NaN also.
		% By Convention norm() returns NaN if one or more of the elements are NaN.
		amaxsum = max(a1, ainf);
		if ( amaxsum >= MyInf ) || isnan(amaxsum)
			b = false;
			PrintResult__('Unstable (sum-norm is too large). ', [], bSuppressUnstable);
			bInfinite = true;
			break;
		end

		% The trace of current A == trace^2 of previous A.
		tr_A = tr_A2;

		% Calculate trace of A*A.
		diag_A2 = sum(A'.*A, 1);
		tr_A2 = full(sum(diag_A2));

		% Get lower bound for spectral radius.
		srlow_A_init = abs(tr_A)/n;

		% For these more refined lower bounds, see:
		%	"Lower bounds for spectral radius of a matrix
		%	- Bill G. Horne, 1995, NECI Technical report 95-14.
		rx = (n*tr_A2 - (tr_A)^2);
		if rx > 0
			srlow_A = srlow_A_init + sqrt(rx/(n-1))/n;
		elseif tr_A2 < 0
			srlow_A = sqrt( (tr_A^2 - tr_A2)/n/(n-1) );
		else
			srlow_A = srlow_A_init;
		end

		% Get spectral radius lower bound for A^2.
		% Saves one expensive A*A multiplication.
		srlow_A2 = abs(tr_A2)/n;

		% Check both bounds. One^2 since we are really checking the next step.
		One = MyOne^loopexp;
		if (srlow_A > One)
			PrintResult__('Unstable (Trace A is too large). ', [], bSuppressUnstable);
			b = false;
			break;
		elseif (srlow_A2 > One^2)
			PrintResult__('Unstable (Trace A-square is too large). ', [], bSuppressUnstable);
			b = false;
			break;
		end

		% Experimental.
		if bUseHeuristic && ( srlow_A_init + (n-1)*(srlow_A - srlow_A_init)/2 > One )
			PrintResult__('Unstable (Heuristic) (Extrapolated bound too large). ', [], bSuppressUnstable);
			b = false;
			bHeuristic = true;
			break;
		end

		loopexp = 2*loopexp;
	end
end

%%
function M = SquareVAR1TypeMatrix__(M, ktop)
	n = length(M);
	M = cat( 1, M(1:ktop, :)*M, M(1 : n-ktop, :) );
end

%%
function M = MultByVAR1TypeMatrix__(M, B, ktop) %#ok<DEFNU>
	n = length(M);
	M = cat( 1, M(1:ktop, :)*B, B(1 : n-ktop, :) );
end

%%
function [bstable, checkfailed, sr] = CheckSpectralRadius__(Op)

	print_x('%s', 'Checking spectral radius ...  ');
	
	[sr, checkfailed] = Op.SpectralRadius();
	if checkfailed
		sr = NaN;
		bstable = false;
		
	else
		bstable = ( sr < 1-eps );
	end
	
end

%%
function PrintResult__(msg, fmt, bSuppressDisplay)
	if (nargin < 3), bSuppressDisplay = false; end
	if (nargin < 2) || isempty(fmt), fmt = '%s\n'; end
	if ~bSuppressDisplay
		print_x(fmt, msg);
	end
end
%%
