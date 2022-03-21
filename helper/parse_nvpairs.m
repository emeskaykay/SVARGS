% [params, remaining] = parse_nvpairs(nvpairs, varargin)
% Parse user input that is in the form of name value pairs and return a structure
% having the fieldnames and values specified in the input.
%
% nvpairs:
%	This is in the form of a cell array
%		{variable_name_1, value_1, variable_name_2, value_2, ... }
%	Typically this is the varargin argument to a function that expects name-value pairs.
%
% varargin [optional]:
% A list of specifications whose elements are cell arrays of the form:
%	{expected_name, expected_type, default_value, [constraint]}
%	expected_name:
%		A string variable name that is expected to be present in nvpairs.
%		** If a name is present in nvpairs but _not_ present in the specification list, an
%			error is thrown unless the second output argument (remaining) is present. If
%			'remaining' is present, then unrecognized names, paired with their values, will
%			be returned in 'remaining'.
%		** Any name present the specification list but not present in nvpairs, will have its
%			corresponding struct field created and assigned the value default_value.
%		** If the specification list is empty or not specified, then all names and values
%			in nvpairs are extracted and assigned empty values.
%
%	expected_type:
%		classname or builtin type corresponding to expected_name.
%		To omit type checking for the particular name, set expected_type to an empty string.
%
%	default_value:
%		Default value to assign to the variable if it is not in the set of name-value pairs.
%		If empty, then the default value is the empty value ([]).
%
%	constraint [optional];
%		Handle to a function or anonymous function or string form of an anonymous function.
%		The function signature should be:
%				[ok, msg] = constraint(x, other).
%		Where "msg" and "other" are optional in the signature.
%
%		The first input to constraint is the value of the variable whose name is <expected_name>.
%		The second (optional) input is a struct whose fieldnames are those names that appear prior to
%		<expected_name> in varargin and which have the corresponding value. The field values in
%		'other' will have been assigned defaults if needed and also will have been constraint checked.
%
%		The constraint function must return a logical value indicating whether the constraint
%		is satisfied. If the constraint is violated, the function can optionally return a second
%		string output indicating the nature of the constraint violation. The first constraint
%		violation will result in an error with the corresponding message printed.
%
%		Simple examples:
%			@is_positive_integer (here is "is_positive_integer" is just a suggestive name)
%			'@(x, ~) x > 0'
%			( @(x, c) (-1 < x) && (x < 5) && (x < c.somename) )
%
% A variable cannot have more than one specification.
%
% RETURNS:
% params:
%	A structure with field names and values.
%
% remaining:
%	If the second output argument (remaining) is present, then unrecognized names (ie:
%	names in nvpairs that were not present in the specification list) along with their
%	values are returned as a cell array of name-value pairs. If the second output argument
%	is NOT present and remaining is non-empty, then an error is thrown.
%
% NOTE: Comparison is case sensitive.
% See also: parse_nvpairs_copy_defaults.m
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
function [params, remaining] = parse_nvpairs(nvpairs, varargin)
	np = numel(nvpairs)/2;
	assert(iscell(nvpairs), 'Input must be a cell array of name value pairs.');
	assert(rem(np,1) == 0, 'Cell array input must be in the form of Name-Value pairs.');
	
	namelist = nvpairs(1:2:end);
	valuelist = nvpairs(2:2:end);
	
	if isempty(varargin)
		assert(iscellstr(namelist) && all(cellfun(@isvarname, namelist)), ...
			'Each pair must consist of a valid variable name and a variable value.');
		
		params = struct;
		for i = 1 : np
			params.(namelist{i}) = valuelist{i};
		end
		remaining = {};
		
		return;
	end
	
	trivial_constraint = (@(~,~) true);
	bfound = false(np, 1);
	expected_namelist = {};
	for i = 1 : numel(varargin)
		spec = varargin{i};
		assert( iscell(spec) && (numel(spec) >= 3) && (numel(spec) <= 4), ['parse_nvpairs: ', ...
			'Each specification must be a cell of the form ', ...
			'{<varname>, <vartype>, <defval>, <constraint>(optional)}.']);

% 		% Get rid of type of default value for the previous iteration to prevent issues
% 		% with inplicit conversion. Example: if default_value was an empty 'double' and
% 		% spec{3} is a 'char', then the assignment default_value = spec{3} will produce
% 		% an empty double instead of an empty char.
% 		clear default_value
		
		expected_name = spec{1};
		expected_type = spec{2};
		default_value = spec{3};
		if numel(spec) < 4
			spec{4} = trivial_constraint;
		end
		
		assert_name_msg__(isvarname(expected_name), expected_name, ' Expected name must be a valid variable name.');
		assert_name_msg__(isempty(expected_type) || isvarname(expected_type), ...
			expected_name, 'Expected type must be a valid type string.');

		assert_name_msg__(isempty(expected_namelist) || ~any(strcmp(expected_namelist, expected_name)), ...
			expected_name, 'Has more than one specification.');
		
		if isvarname(expected_type) && ~isempty(default_value)
			assert(isa(default_value, expected_type), ...
				expected_name, ['Type of default value did not match ', expected_type]);
		end
		
		fconstraint = spec{4};
		if ischar(fconstraint), fconstraint = str2func(fconstraint); end
		assert_name_msg__(isa(fconstraint, 'function_handle'), ...
			expected_name, 'Constraint must be a function string or function handle');
		
		b = vecc(strcmp(namelist, expected_name));
		r = nnz(b);
		if r == 0
			params.(expected_name) = default_value;
			
		elseif r > 1
			assert_name_msg__(false, expected_name, 'Name-value pairs cannot be repeated.');
			
		else % r == 1
			assert_name_msg__(isempty(expected_type) || isa(valuelist{b}, expected_type), ...
				namelist{b}, ['Type did not match ', expected_type]);
			params.(expected_name) = valuelist{b};
			bfound = (bfound | b);
		end
		
		fcstr = func2str(fconstraint);
		
		nfin = nargin(fconstraint);
		assert_name_msg__(nfin >= 1, expected_name, 'Constraint must take at least one input.');
		
		if (fcstr(1) == '@')
			nfout = 1;
		else
			nfout = nargout(fconstraint);
		end
		assert_name_msg__(nfout >= 1, expected_name, 'Constraint must return at least one output.');
		
		if (nfout == 1)
			if nfin == 1
				ok = fconstraint(params.(expected_name));
				
			else
				ok = fconstraint(params.(expected_name), params);
			end
			
		else
			if nfin == 1
				[ok, cmsg] = fconstraint(params.(expected_name));
				
			else
				[ok, cmsg] = fconstraint(params.(expected_name), params);
			end
			
			disp(cmsg);
		end
		
		if false == ok
			error('parse_nvpairs:constraint', [expected_name, ': Constraint was violated.']);
		end
		
		expected_namelist = [expected_namelist, {expected_name}]; %#ok<AGROW>
	end
	
	bunknown = ~bfound;
	remaining = vecr( [namelist(bunknown); valuelist(bunknown)] );
	if (nargout == 1) && ~isempty(remaining)
		error('parse_nvpairs: Unrecognized names found. ');
	end
end
%%

function msg = assert_name_msg__(b, name, msg)
	msg = ['parse_nvpairs: ', name, ': ', msg];
	assert(b, msg);
end

%%
