%
%
%
%
%
%% Handle class
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
classdef Svargs_hist < dynamicprops
	
	properties (Access = 'private')
		eqnstathist_vec;
		hcriterion;
	end
	
	properties (SetAccess = 'private')
		ntimepts_used;
	end
	
	methods
		function This = Svargs_hist(ntim, hcrit)
			This.ntimepts_used = ntim;
			This.hcriterion = hcrit;
		end
		
		function eqsh = NewEquationHist(This, C0, Sig0)
			kcurr = numel(This.eqnstathist_vec);
			eqsh = Svargs_hist_eqn(kcurr+1, C0, Sig0, This.ntimepts_used, This.hcriterion);
			This.eqnstathist_vec = [This.eqnstathist_vec, eqsh];
		end
		
		% Accesssors.
		function r = EquationCount(This)
			r = numel(This.eqnstathist_vec);
		end
		
		function sel = Selector(This, iblk)
			sel = This.GetEquationHist(iblk).Selector;
		end
		
		function eqsh = GetEquationHist(This, iblk)
			if nargin < 2
				eqsh = This.eqnstathist_vec(end);
			else
				eqsh = This.eqnstathist_vec(iblk);
			end
		end
		
		function [Coeffs, CoeffsInitial] = GetBestCoeffs(This, startcol, endcol)
			neqblocks = numel(This.eqnstathist_vec);
			assert(neqblocks > 0, 'Svargs_hist:GetBestCoeffs, Equation history is empty. ');
			
			if (nargin < 3) || isempty(endcol), endcol = This.eqnstathist_vec(1).ncols; end
			if (nargin < 2) || isempty(startcol), startcol = 1; end
			
			k = This.eqnstathist_vec(1).nrows;
			for iblk = 2 : neqblocks
				assert(This.eqnstathist_vec(iblk).ncols == This.eqnstathist_vec(1).ncols);
				k = k + This.eqnstathist_vec(iblk).nrows;
			end
			
			ncols = endcol - startcol + 1;
			Coeffs(k, ncols) = 0;
			CoeffsInitial(k, ncols) = 0;
			last = 0;
			for iblk = 1 : neqblocks
				blklen = This.eqnstathist_vec(iblk).nrows;
				coeffs = This.eqnstathist_vec(iblk).GetBestCoeffs();
				Coeffs(last+1:last+blklen, :) = coeffs(:, startcol:endcol);
				CoeffsInitial(last+1:last+blklen, :) = This.eqnstathist_vec(iblk).coeffs_initial(:, startcol:endcol);
				last = last + blklen;
			end
		end
		
	end

end
%%
