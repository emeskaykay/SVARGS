% [idbest valbest] = Best(This, hcrit, ids, mintype, varargin)
% Returns the best model id, globally or locally, as per the specified option
% and criterion.
%	
% hcrit [required for optname = 'quality' and for optval = 'best' or 'local']:
%	Handle to criterion method. One of:
%		@aic
%		@aic_c
%		@sbc
%		@sbc_x
%		@hqc
%		@fpe
%	or handle to a method having the same signature as the above.
%
% ids [default: all ids]:
%	A vector of ids from which the best is to be chosen.
%
% mintype [default: 'global']:
% This must be one of:
%	'global': Return model id and quality value of the global best model.
%
%	'local': Return model id and quality value of best model which is also
%	a turning point in the quality vs. complexity graph. If there is no such model,
%	then the returned id is empty.
%
% varargin:
% Name-Value  pairs
% Name: reltol, Value: numeric tolerance for quality [default = 0]:
%	Relative criterion tolerance. Relevant only when optname = 'best'. Quality
%	values lying in the range [minval, minval+rangewidth*reltol)] are considered
%	for determining the best id(s). Here minval is the best quality value, and
%	rangewidth is size of the range of quality values.
%
% Name: complexity, Value: 'lowest', 'highest', 'average' or 'all'
%	In the case where there is more than one best value value (valbest say), this
%	specifies which id to return. Specifying 'lowest' returns the id with the lowest
%	complexity, 'highest' returns the id with the highest complexity and 'average'
%	returns the id closest to the average complexity of the ids whose quality equals
%	valbest. Specifying 'all' returns all the ids having quality valbest.
%
% RETURNS:
% idbest:
%	Best id as per options.
% 
% valbest:
%	The criterion value (quality) corresponding to the id.
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
function [idbest, valbest] = Best(This, hcrit, ids, mintype, varargin)
	
	opts = parse_nvpairs(varargin, {'reltol', 'double', 0}, {'complexity', 'char', 'lowest'});
	
	if (nargin < 4) || isempty(mintype), mintype = 'global'; end
	if (nargin < 3) || isempty(ids), ids = This.idrecmap.keys(); end
	
	reltol = opts.reltol;
	cxpos = opts.complexity;
	
	assert(isscalar(reltol) & (reltol >= 0));
	assert(ischar(mintype));
	
	if iscell(ids)
		ids = [ids{:}];
	end
	
	assert(isa(hcrit, 'function_handle'));
	assert(isnumeric(ids) && isvector(ids));

	if strcmpi(mintype, 'global')
		vals = Quality(This, hcrit, ids);
		maxval = max(vals);
		[idbest, valbest] = findbestids__(This, vals, ids, cxpos, reltol, maxval);

	elseif strcmpi(mintype, 'local')
		[~, ic] = sort(Complexity(This, ids), 'ascend');
		ids = ids(ic);
		ids = vecc(ids);
		
		allvals = Quality(This, hcrit, ids);
		allvals = vecc(allvals);

		bgood = ~(isnan(allvals) | isinf(allvals));
		allvals = allvals(bgood);
		ids = ids(bgood);
		
		idbest = [];
		valbest = [];
		[~, bmin] = extrema(allvals, 1);
		if any(bmin)
			vals = allvals(bmin);
			ids = ids(bmin);
			maxval = max(allvals);
			[idbest, valbest] = findbestids__(This, vals, ids, cxpos, reltol, maxval);
		end
		
	else
		assert(false, ['Unknown option: ', mintype]);
	end
	
end
%%

function [idbest, valbest] = findbestids__(This, vals, ids, cxpos, reltol, maxval)
	valmin = nanmin(vals);
	width = maxval - valmin;
			
	% Should work correctly even if there is only one id.
	ids = ids((vals <= (valmin + width*reltol)));

	cxvec = Complexity(This, ids);
	switch (lower(cxpos))
		case 'lowest'
			[~, i] = nanmin(cxvec);
			idbest = ids(i);
			valbest = valmin;
			
		case 'highest'
			[~, i] = nanmax(cxvec);
			idbest = ids(i);
			valbest = valmin;
			
		case 'average'
			cxavg = nanmean(cxvec);
			[~, i] = min(abs(cx - cxavg));
			valbest = ids(i);
			
		case 'all'
			[~, ixperm] = sort(cxvec);
			valbest = vals(vals(ixperm) == valmin);

		otherwise
			assert(false);
	end
end

%%
