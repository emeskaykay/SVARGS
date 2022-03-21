% print_obj(fileid, object)
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
function print_obj(fileid, object)
	if (nargin == 1), object = ''; end
	if (nargin == 0), object = ''; fileid = 1; end
	
	if ischar(object)
		outstr = sprintf('%s', object);
		if length(outstr) > 100
			outstr = outstr(1:97);
			outstr(98:100) = '...';
			outstr(101:end) = [];
		end
		fprintf(fileid, outstr);
		fprintf(fileid, '\n');
		
	else
		if isnumeric(object)
			outstr = evalc('disp(object)');
		else
			outstr = evalc('display(object)');
		end
		fprintf(fileid, outstr);
	end

end
%%
