% denude_coeffs(coeffs_dense, numcoeffs, csizemin, srlims, buseheuristic, bfastabort)
%
% coeffs_dense:
%	(k x k x p) dense coefficient array.
%
% numcoeffs:
%	Approximate number of coefficients required. The returned sparse
%	LagOperator will have _approximately_ this number of coefficients.
%
% csizelims:
%	Pair of values, [csizemin, csizemax], specifying minimum and maximum absolute coefficient size.
%
% srmax:
%	Maximum value of spectral radius.
%
% Returns:
% opsparse:
%	Sparse LagOperator object.
%	
% bstable:
%	Logical value that is true if operator corresponding to C0 is stable.
%
% sr:
%	Spectral radius of operator.
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
function [opsparse, found, sr, niter] = denude_coeffs(coeffs_dense, numcoeffs, csizelims, srmax, bfastabort, verbose)
	[k, ~, p] = size(coeffs_dense);
	
	dlags = 1:p;
	cmax = max(abs(coeffs_dense(:)));
	coeffcnt = nnz(coeffs_dense);
	
	if (nargin < 5) || isempty(bfastabort), bfastabort = false; end
	if (nargin < 4) || isempty(srmax), srmax = 1 - (1+rand)/k/p; end
	if (nargin < 3) || isempty(csizelims), csizelims = [0, Inf]; end
	
	assert( isnumeric(csizelims) && (numel(csizelims) == 2) );
	csizemin = csizelims(1);
	csizemax = csizelims(2);
	assert((csizemin >= 0) && (csizemax >= 0.5), 'VARSparseViaDense: Invalid coefficient size limits. ');
	assert((csizemin < 0.5), 'VARSparseViaDense: Minimum coefficient size is too large for efficient search. ');
	
	% Try using true only if stable model is not found.
	tapered_deflate = false;
	deflator = 0.99;
	
	found = false;
	elapsed_total = 0;
	sr = Inf;
	niter = 0;
	abort = false;
	mincmax = csizemin + 0.5;
	
	if false == verbose
		dispmasker = Display_masker('char');
		olddisp = dispmasker.Flag('char', true);
	end
	
	progress = Progress_display('Denuding coeffs', {'Iter', 'Max_size', 'SR', 'C_cnt', 'Ord'});
	while ( ~abort && ~found )
		tic;

		C0 = zero_out__(coeffs_dense, dlags, csizemin);
		if tapered_deflate
			C0 = deflate__(C0, dlags, deflator^niter, csizemin);
		else
			C0 = deflate__(C0, 1, deflator^niter, csizemin);
		end
		
		% Denude coefficients uniformly to get numcoeffs coeffs.
		ixvec = find(C0 ~= 0);
		ix = datasample(ixvec, numcoeffs, 'replace', false);
		B = true(size(C0));
		B(ix) = false;
		C0(B) = 0;
		
		cmax = max(abs(C0(:)));
		if (cmax < csizemin+0.5)
			break; end
		
		opsparse = LagOperator(dlags, C0);
		[found, ~,~, loopexp] = opsparse.IsStable();
		
		ord = opsparse.Order();
		
		coeffcnt = nnz(C0);
		progress.update([niter, cmax, sr, coeffcnt, ord]);
		
		if ~found && (loopexp >= 2^(ord-1))
			sr = opsparse.SpectralRadius();
			progress.update([niter, cmax, sr, coeffcnt, ord]);
			r = srmax/sr;
			if (r > 1)
				disp('Either SpectralRadius() computation was inaccurate or IsStable returned false negative.');
				disp('Redoing loop.');
				continue;
			end
			t = r;
			while ~found && (t^ord > 0.8)
				C0 = deflate__(C0, dlags, r, csizemin); % Taper to guarantee required reduction in sr.
				opsparse = LagOperator(dlags, C0);
				found = opsparse.IsStable();
				cmax = max(abs(C0(:)));
				ord = opsparse.Order();
				progress.update([niter, cmax, sr, coeffcnt, ord]);
				r = min(srmax/sr, r);
				t = t*r;
			end
		end
		
		found = found && (cmax >= mincmax);
		
		elapsed = toc;
		elapsed_total = elapsed_total + elapsed;
		abort = ~found && bfastabort && ((elapsed > 5) || (elapsed > 3*(elapsed_total/niter)));
		niter = niter + 1;
	end
	
	if found
		if cmax > csizemax
			r = csizemax/cmax;
			C0 = deflate__(C0, 1, r, csizemin); % No need to taper.
			opsparse = LagOperator(dlags, C0);
		end
		sr = opsparse.SpectralRadius();
		progress.complete([niter, cmax, sr, coeffcnt, ord]);
	else
		progress.complete([niter, cmax, sr, coeffcnt, ord]);
		disp('Not found. ');
	end
	
	if false == verbose
		dispmasker.Flag('char', olddisp);
	end
end

%%
function C = shift_out__(C, csizemin)
	B = (C > -csizemin) & (C < 0);
	C(B) = C(B) - (csizemin + max(C(B)));

	B = (C > 0) & (C < csizemin);
	C(B) = C(B) + (csizemin - min(C(B)));
end

%%
function C = zero_out__(C, dlags, csizemin)
	margin = 0.02;
	taper = (1-margin).^dlags;
	csizemin = csizemin./taper;
	csizemin = reshape(csizemin, [1,1,numel(dlags)]);
	B = bsxfun(@gt, C, -csizemin) & bsxfun(@lt, C, csizemin);
	C(B) = 0;
end

function C0 = deflate__(C, dlags, r, csizemin)
	C0 = C;
	if r < 1
		r_vec = reshape(r.^dlags, 1,1,[]);
		Cplus = csizemin + bsxfun(@times, (C - csizemin), r_vec);
		Cminus = -csizemin + bsxfun(@times, (C + csizemin), r_vec);
		C0(C > csizemin) = Cplus(C > csizemin);
		C0(C < -csizemin) = Cminus(C < -csizemin);
	end
end

%%
