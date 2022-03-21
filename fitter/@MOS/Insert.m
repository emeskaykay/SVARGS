% Insert
% Records the model complexities and costs corresponding to the specified model id(s).
%
% ids:
%	Vector of numeric model ids to insert.
%
% complexities:
%	Vector of numeric model complexities.
%
% errcovs:
%	Array of error covariance matrices. 3-dimensional array of slices.
%	Only the log absolute determinant of the covariance matrix is recorded.
%
% owwarn:
%	If true, warns if any id's were overwritten.
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
function Insert(This, ids, complexities, errcovs, owwarn)

	if nargin < 5 || isempty(owwarn), owwarn = true; end
	idcount = numel(ids);
	
	assert( isnumeric(ids) && isvector(ids) );
	assert( isnumeric(complexities) && isvector(complexities) );
	assert( isnumeric(errcovs) );
	assert( (size(errcovs,1) == size(errcovs,2)), 'MOS:Insert: ''errcovs'' must be a 3D numeric array. ' );
	assert( numel(complexities) == idcount );
	
	overwritten = 0;
	errcovs = num2cell(errcovs, [1,2]);
	for j = 1 : idcount
		overwritten = overwritten + This.idrecmap.isKey(ids(j));
		This.idrecmap(ids(j)) = MOS_record(complexities(j), errcovs{j});
	end
	
	if owwarn && (overwritten > 0)
		disp( ['MOS:Insert : ', num2str(overwritten), ' record(s) overwritten. '] );
	end

end
%%
