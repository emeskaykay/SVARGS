% qvals = Quality(This, hCriterion, ids)
% For each of the specified id's, returns the model quality value (ic value)
% as per the specified information criterion (ic).
%
% hCriterion:
% One of:
%	@aic
%	@aic_c
%	@sbc
%	@sbc_x
%	@hqc
%	@fpe
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
function qvals = Quality(This, hCriterion, ids)

	assert( isnumeric(ids) && isvector(ids) );
	
	Recs = This.idrecmap.values(num2cell(ids));
	Recs = [Recs{:}];
	
	costs = reshape([Recs.cost], size(ids));
	complexities = reshape([Recs.complexity], size(ids));
	qvals = arrayfun(hCriterion, repmat(This, size(ids)), costs, complexities, 'UniformOutput', true);

end
%%
