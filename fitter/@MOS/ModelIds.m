% [ids vals] = ModelIds(This, sortby, hcrit)
% Returns a vector of model id's as per the specified sort order and corresponding
% values. ModelIds() without any arguments returns a vector of model id's sorted by id.
%
% sortby:
%	This this must be one of:
%		'sort_quality'		: Return id's and quality values sorted by model quality.
%								In this case hcrit must be specified.
%		'sort_complexity'	: Return id's and quality values sorted by model complexity.
%		'sort_cost'			: Return id's and costs sorted by model cost.
%		'sort_modelid'		: Return id's sorted by id. (Same as no-argument call).
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
% RETURNS:
% ids:
%	Vector of ids as per sort option.
% 
% vals:
%	The values of the of the type specified in the sortby argument.
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
function [ids, vals] = ModelIds(This, sortby, hcrit)
	
	ids = This.idrecmap.keys();
	ids = [ids{:}];
	
	if (nargin <= 2) || strcmpi(sortby, 'modelid')
		vals = [];
		return;
	end
	
	if strcmpi(sortby, 'quality')
		assert( (nargin == 3) && isa(hcrit, 'function_handle') );
		vals = Quality(This, hcrit, ids);

	elseif strcmpi(sortby, 'complexity')
		vals = Complexity(This, ids);

	elseif strcmpi(sortby, 'cost')
		vals = Cost(This, ids);
	end

	[vals, ix] = sort(vals, 'ascend');
	ids = ids(ix);
	
end
%%
