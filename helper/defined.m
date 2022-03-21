% defined
% Checks if a variable with specified name exists in the specified workspace
% and takes on any one of the specified values. If no value is specified, the
% function simply checks if the variable exists (even if empty) in the specified
% workspace.
%
% Usage: defined('...')
% For example: if defined('debugging')
%					checks if there is a variable called 'debugging in the base workspace.
%
%			or  : if defined('mode', [1 3 5])
%					checks if a variable called 'mode' is defined and equal to 1 or 3 or 5.
%
%			or  : if defined({'mode', debug', 'suppress'}, {[1,3,5], [true, false], false});
%					checks if a variable called 'mode' is defined and equal to 1 or 3 or 5.
%					checks if a variable called 'debug' is defined and is either true or false.
%					checks if a variable called 'suppress is defined and has value false.
%
%			or  : if defined({'mode', debug', 'level', 'suppress'}, {[1,2,5], [true, false], [], false});
%					Same as above and also checks if 'level' is defined. The [] signifies that the value
%					should not be checked.
%
% varstrings:
%	Either a string representing the name of variable to check
%	OR a cell array of strings each of which is to be checked.
%
% wkspace [default: 'base workspace ('base')]
%	The workspace where the variable is to be found.
%
% values [optional]
%	Either a scalar/vector of values OR a cell array of scalar/vector values
%	corresponding to each name in varstrings.
%
% Returns:
% b:
%	Boolean value indicating whether the base workspace variable specified
%	exists and is equal to one of the specified values.
%
% bvec:
%	Boolean vector indicating whether each base workspace variable specified
%	exists and is equal to one of the specified values.
%
% varvals:
%	The actual value of the specified variable. If a cell array of names was
%	provided in the varstring argument, then this is a cell array of the same
%	length containing the actual values of each of the variables specified.
%	Non-existent variables will result in an empty entry in varvals.
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
function [b, bvec, varvals] = defined(varstrings, wkspace, values)

	if (nargin < 3), values = {}; end
	if (nargin < 2) || isempty(wkspace), wkspace = 'base'; end
	
	cellform = true;
	assert( ~isempty(varstrings) );
	if ischar(varstrings)
		cellform = false;
		varstrings = {varstrings};
	end
	if ~iscell(values), values = {values}; end
	
	varvals = cell(numel(varstrings), 1);
	
	bvec = false(numel(varstrings), 1);
	for i = 1 : numel(bvec)
		assert(ischar(varstrings{i}), 'defined: First argument must be a string or cell array of strings. ');
		if evalin(wkspace, ['exist(', '''', varstrings{i}, '''', ', ''var'')'])
			bvec(i) = true;
			varvals{i} = evalin(wkspace, varstrings{i});
			if ~isempty(values)
				bvec(i) = bvec(i) && (isempty(values{i}) || any(values{i} == varvals{i}));
			end
		end
	end
	b = all(bvec);
	
	if ~cellform
		varvals = varvals{1};
	end
end
%%
