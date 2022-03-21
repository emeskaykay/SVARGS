% Svargs_hist_eqn
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
classdef Svargs_hist_eqn < handle
	
	properties (SetAccess = 'private')
		%% Constants.
		objid;
		datalen;				% Number of time points in data.
		nrows;				% Number of rows in coefficient matrix.
		ncols;				% number of columns in coefficient matrix
		coeffs_initial;	% Starting coefficients (default: all zeros).
		sigma_initial;
		
		%% Variables
		% Step number. Changed when new coefficients are entered.
		serial;
		
		% Current level number. Changes are recorded at a given level, overwriting previous
		% changes. When the level is incremented, the model quality is recorded and a new
		% level created. A given level may consist of more than one serial.
		level;
		
		% MOS (Model Order Selector) object. % If required, this can be used later
		% with a different information criterion than the one supplied in the constructor.
		% Modified when level is incremented.
		Selector;
	end
	
	properties (Access = 'private')
		leveldata;	% struct array.
		
		% For tracking the minima of the ic values.
		hcriterion;
		level_local_best;
		quality_local_best;
		level_global_best;
		quality_global_best;
	end
	
	properties (Dependent, SetAccess = 'private')
		coeff_count;	% Coeff count at current level.
		
		% Fields stored in each leveldata struct.
		coeffs;		% Coefficients at current level and serial.
		entryid;		% Id specified by the caller at each entry. This can be of any user defined type.
		sigma;		% Noise covariance matrix at current level and serial.
		cxdelta;		% cxdelta is the change in number of non-zero elements at the last step.
	end
	
	methods (Access = 'private')
		function SetCurrentLevel__(This, lev)
			assert(numel(This.leveldata) >= lev);
			assert((lev == 0) || ~isempty(This.leveldata(lev)));
			This.level = lev;
		end
		
		% Increment level and initialize values to values at previous level.
		function IncrementLevel__(This)
			This.Selector.Insert(This.level, This.coeff_count, This.sigma);
			
			quality_curr_lev = This.Selector.Quality(This.hcriterion, This.level);
			
			if isempty(This.level_local_best) || (quality_curr_lev < This.quality_local_best)
				[This.level_local_best, This.quality_local_best] = This.Selector.Best(This.hcriterion, [], 'local');
			end
			
			if isempty(This.level_global_best) || (quality_curr_lev < This.quality_global_best)
				This.level_global_best = This.level;
				This.quality_global_best = quality_curr_lev;
			end

			lev = This.level + 1;
			
			This.leveldata(lev).entryid = This.entryid;
			This.leveldata(lev).coeffs = This.coeffs;
			This.leveldata(lev).sigma = This.sigma;
			This.leveldata(lev).cxdelta = This.cxdelta;
			
			SetCurrentLevel__(This, lev);
		end
	end
	
	methods
		function v = get.coeff_count(This)
			v = nnz(This.coeffs);
		end
		
		% Access fields stored in leveldata struct.
		function v = get.entryid(This)
			if This.level > 0
				v = This.leveldata(This.level).entryid;
			else
				v = [];
			end
		end
		
		function v = get.coeffs(This)
			if This.level > 0
				v = This.leveldata(This.level).coeffs;
			else
				v = This.coeffs_initial;
			end
		end
		
		function v = get.sigma(This)
			if This.level > 0
				v = This.leveldata(This.level).sigma;
			else
				v = This.sigma_initial;
			end
		end
		
		function v = get.cxdelta(This)
			if This.level > 0
				v = This.leveldata(This.level).cxdelta;
			else
				v = 0;
			end
		end
	end
	
	methods
		function This = Svargs_hist_eqn(iobj, c_init, sig_init, N, hcrit)
			assert(isnumeric(c_init) && ismatrix(c_init));
			assert(isnumeric(sig_init) && ismatrix(sig_init));
			
			[nr, nc] = size(c_init);
			assert( all(size(sig_init) == [nr,nr]) );
			
			This.objid = iobj;
			This.datalen = N;
			This.hcriterion = hcrit;
			This.nrows = nr;
			This.ncols = nc;
			This.coeffs_initial = sparse(c_init);
			This.sigma_initial = sig_init;
			
			This.level = 0;
			This.serial = 0;
			This.Selector = MOS(N, 1, nc); % Note that MOS ia a handle class.
			
			This.IncrementLevel__();
		end
		
		function [levs, quals] = IcTail(This)
			levs = [];
			quals = [];
			last_lev = This.level - 1;
			if ~isempty(This.level_global_best) && (This.level_global_best < last_lev)
				levs = (This.level_global_best+1 : last_lev);
				quals = This.Selector.Quality(This.hcriterion, levs);
			end
		end
		
		function m = IcTailLength(This)
			m = 0;
			if ~isempty(This.level_global_best)
				last_lev = This.level - 1;
				m = last_lev - This.level_global_best;
			end
		end
		
		% Get the best level and quality (as per criterion specified in constructor).
		function [lev, quality] = GetBestLevel(This)
			lev = This.level_global_best;
			quality = This.quality_global_best;
		end
		
		% Get the best coefficient set (as per criterion specified in constructor).
		function [cc, sig, lev, quality] = GetBestCoeffs(This)
			[lev, quality] = This.GetBestLevel();
			
			if isempty(lev)
				cc = [];
				sig = [];
				
			elseif lev > 0
				cc = This.leveldata(lev).coeffs;
				sig = This.leveldata(lev).sigma;
				
			else % lev == 0
				cc = This.coeffs_initial;
				sig = This.sigma_initial;
			end
		end
		
		function c = GetCoeffs(This, lev)
			if lev == 0
				c = This.coeffs_initial;
			else
				assert(lev > 0);
				c = This.leveldata(lev).coeffs;
			end
		end
		
		% Add new coefficients and update existing coefficients.
		% id: id for the entry.
		%
		% ixcols_new, cc_new: (column) indices whose coefficients are to be changed
		% from 0 or set to 0 and the corresponding coefficients. Non-empty arguments
		% here will affect the location indicator.
		%
		% ixcols_updt, cc_updt: (column) indices that are to be modified (existing
		% coeffs at these locations must be non-zero) and the corresponding updated
		% coefficients (which must also be non-zero). These arguments dont affect
		% the location indicator.
		% 
		% sto: noise for the "to" equation.
		function EnterCoeffs(This, id, ixcols_new, cc_new, ixcols_updt, cc_updt, sto)
			lev = This.level;
			This.serial = This.serial + 1;
			
			% Enter new coeffs or zero-out existing
			cc0 = This.leveldata(lev).coeffs(:, ixcols_new);
			assert( all(cc0(:) == 0) || all(cc_new(:) == 0) );
			
			This.leveldata(lev).coeffs(:, ixcols_new) = cc_new(:,:);
			This.leveldata(lev).entryid = id;
			This.leveldata(lev).sigma = sto;
			This.leveldata(lev).cxdelta = nnz(cc_new) - nnz(cc0);
			
			% update existing coeffs.
			cc0 = This.leveldata(lev).coeffs(:, ixcols_updt);
			%assert( all( vecc(cc_updt ~= 0) == vecc(cc0 ~= 0) ) );
			This.leveldata(lev).coeffs(:, ixcols_updt) = cc_updt;
			
			This.IncrementLevel__(); % Record the coeffs and quality.
		end
		
	end
	
end
%%
