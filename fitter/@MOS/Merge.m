% Merge
% Merge two selector objects.
% Usage: sleft.Merge(sright)
%				OR
%			Merge(sleft, sright)
% 
% For any common model ids present, the covariances in 'sleft' object
% are overwritten by the covariances in 'sright'. The object 'sright' is
% not modified.
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
function Merge(sleft, sright)

	idsright = sright.idrecmap.keys();
	b = sleft.idrecmap.isKey(idsright);
	
	% If duplicates are not removed the stupid thing excludes both values.
	sleft.idrecmap.remove(idsright(b));
	
	irmap1 = sleft.idrecmap;
	irmap2 = sright.idrecmap;
	sleft.idrecmap = [irmap1; irmap2]; % Notation for merging maps.
	
end
%%
