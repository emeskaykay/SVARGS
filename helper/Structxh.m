% Structxh handle class
% struct with extra features and locking.
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
classdef Structxh < handle
	properties (Access = private)
		sdata;
		islocked = false;
	end
	
	methods (Static, Access = private)
		function showstruct__(prefix, s)
			if isempty(s)
				return;
			end
			
			fnames = fieldnames(s);
			if isempty(fnames)
				disp('Structxh: Empty Structxh object. ');
				return;
			end
			
			for i = 1 : numel(fnames)
				val = s.(fnames{i});
				if ~isempty(prefix)
					prefix2 = [prefix, '.', fnames{i}];
				else
					prefix2 = fnames{i};
				end
				
				if isstruct(val)
					Structxh.showstruct__(prefix2, val);
					
				else
					disp([prefix2, ': ']);
					if ischar(val)
						disp(['    ', val, newline]);
						
					elseif iscell(val)
						for v = vecr(val)
							disp(['    ', v{1}]);
						end
						disp(newline);
						
					else
						disp(vecc(val));
					end
				end
			end
		end
		
		% Descend to numeric value or object.
		function [leaf, par, fnames] = get_leaf__(s, dotpath)
			assert(isstruct(s) && isscalar(s));
			assert(ischar(dotpath) && ~isempty(dotpath));

			fnames = strsplit(dotpath, '.');
			
			leaf = s;
			depth = numel(fnames);
			i = 0;
			while (i < depth)
				par = leaf;
				if isstruct(leaf) && isfield(leaf, fnames{i+1})
					i = i + 1;
					leaf = leaf.(fnames{i});
				else
					break;
				end
			end
		end

		function n = count_fields__(s, depth)
			assert(depth >= 0);
			
			n = 1;
			if (depth > 0) && isstruct(s)
				fnames = fieldnames(s);
				for i = 1 : numel(fnames)
					n = n + Structxh.count_fields__(s.(fnames{i}), depth-1);
				end
			end
		end
		
		function n = count_leaf_values__(s, depth, omit_empty_values)
			assert((omit_empty_values == 1) || (omit_empty_values == 0))
			assert(depth >= 0);
			
			% if s is not a struct it is assumed to be a value (possibly empty).
			if ~isstruct(s)
				n = ~(isempty(s) && omit_empty_values);
			
			elseif depth == 0
				n = 0;
				
			else
				n = 0;
				fnames = fieldnames(s);
				for i = 1 : numel(fnames)
					n = n + Structxh.count_leaf_values__(s.(fnames{i}), depth-1, omit_empty_values);
				end
			end
		end
	end
	
	methods (Access = private)
		function check_locked__(hndl, dotnames)
			ok = false;
			if hndl.islocked
				if (nargin == 2) && ~isempty(dotnames)
					ok = has_field(hndl, dotnames);
				end
				if ~ok
					error('Structxh:struct_locked', 'Structxh object must be unlocked first. ');
				end
			end
		end
		
		function deletef__(hndl, dotnames)
			if ~isfield(hndl.sdata, dotnames{1})
				% do nothing.
				
			elseif numel(dotnames) == 1
				hndl.sdata = rmfield(hndl.sdata, dotnames{1});
				
			else
				h = Structxh(hndl.sdata.(dotnames{1}));
				h.deletef__(dotnames(2:end));
				hndl.sdata.(dotnames{1}) = h.sdata;
			end
		end
		
	end
	
	methods
		function hndl = Structxh(s)
			if (nargin == 0) || isempty(s), s = struct; end
			set_data(hndl, s);
		end
		
		function load_data(hndl, filepath)
			hndl.check_locked__();
			hndl.sdata = load(filepath);
		end
		
		% Get a deep copy.
		function structxhobj = get_copy(hndl)
			structxhobj = Structxh;
			structxhobj.sdata = hndl.sdata;
		end
		
		% Saves all fields except those beginning with 'TEMP.
		function save_data(hndl, filepath)
			s = hndl.sdata;
			if isfield(s, 'TEMP')
				s = rmfield(s, 'TEMP');
			end
			save(filepath, '-struct', 's');
		end
		
		function n = field_count(hndl, depth)
			n = hndl.count_fields__(hndl.sdata, depth);
			n = n - 1; % omit head.
		end
		
		function n = value_count(hndl, depth, omit_empty_values)
			n = hndl.count_leaf_values__(hndl.sdata, depth, omit_empty_values);
		end
		
		function yes = is_locked(hndl)
			yes = hndl.islocked;
		end
		
		function show(hndl, fdotpath)
			if exist('fdotpath', 'var') && ~isempty(fdotpath)
				s = get_field(hndl, fdotpath);
				if isempty(s)
					disp('Structxh: Dot path does not exist. ');
					return;
				end
				
			else
				fdotpath = '';
				s = get_data(hndl);
			end
			Structxh.showstruct__(fdotpath, s);
		end
		
		% When locked no fields can be added via set_field() or removed via delete_field().
		% Also the data cannot be set via set_data() or load_data(). An attempt to add
		% or remove fields or set the data will result in an exception. However, fields
		% can modified. This can also be useful to prevent new fields from being mistakenly
		% added due to typos.
		function lock(hndl)
			hndl.islocked = true;
		end
		
		% Unlocks the handle.
		function unlock(hndl)
			hndl.islocked = false;
		end
		
		% matlab isfield() does not do nested structs.
		% However getfield() and setfield() do handle nested structs. See the docs.
		function [fexists, val] = has_field(hndl, dotnames)
			if ~isfield(hndl.sdata, dotnames{1})
				fexists = false;
				val = [];
				
			elseif numel(dotnames) == 1
				fexists = true;
				if nargout > 1
					val = hndl.sdata.(dotnames{1});
				end
				
			elseif ~isstruct(hndl.sdata.(dotnames{1})) % => spurious dotnames after 1.
				fexists = false;
				val = [];
				
			else
				h = Structxh(hndl.sdata.(dotnames{1}));
				if nargout > 1
					[fexists, val] = has_field(h, dotnames(2:end));
					
				else
					fexists = has_field(h, dotnames(2:end));
				end
			end
		end
		
		% Returns copy of internal data structure.
		function s = get_data(hndl)
			s = hndl.sdata;
		end
		
		% Copies supplied structure to internal data.
		% s must be a struct or a Structxh object.
		function set_data(hndl, s)
			hndl.check_locked__();
			
			if isstruct(s)
				hndl.sdata = s;

			elseif isa(s, 'Structxh')
				hndl.sdata = s.sdata;

			else
				assert(false, 'Structxh:set_data(), input must be a struct or a Structxh.');
			end
		end
		
		% Sets the internal data structure to the empty struct, ie: removes all fields.
		function clear_data(hndl)
			hndl.sdata = struct;
		end
		
		% Checks for (possibly nested) fieldname and returns value.
		% The full (unique) path to the fieldname must be supplied.
		% fdotpath:
		%	fieldname or dotted path consisting of fieldnames
		%	(example 'namex' or 'namex.namey.namez')
		%
		% OUTPUT:
		% val:
		%	value of field data or empty if field does not exist.
		%
		% fexists:
		%	true if field exists.
		function [val, fexists] = get_field(hndl, fdotpath)
			dotnames = strsplit(fdotpath, '.');
			[fexists, val] = has_field(hndl, dotnames);
		end
		
		% Sets value of the field specified by fdotpath to the value
		% specified by val. If the field does not exist, it is created.
		%
		% fdotpath:
		%	fieldname or dotted path consisting of fieldnames
		%	(example 'namex' or 'namex.namey.namez')
		%
		% val:
		%	value to which the field is to be set.
		function set_field(hndl, fdotpath, val)
			dotnames = strsplit(fdotpath, '.');
			hndl.check_locked__(dotnames);
			
			[leaf, ~, fnames] = Structxh.get_leaf__(hndl.sdata, fdotpath);
			if isstruct(leaf)
				hndl.sdata = setfield(hndl.sdata, dotnames{:}, val);
				
			else
				if isobject(leaf)
					disp(['Structxh::set_field, Overwriting user defined object ' fnames{end}]);
				end
				
				n = numel(fnames);
				fnames(end) = [];
				S = struct;
				S = setfield(S, dotnames{n:end}, val); %#ok<SFLD>
				if isempty(fnames)
					hndl.sdata = S;
				else
					hndl.sdata = setfield(hndl.sdata, fnames{:}, S);
				end
			end
		end
		
		% Sets field corresponding to fdotpath to empty value ([]).
		% If the field does not exist, it is created and set to empty.
		% If the field does not exist, the handle must be unlocked.
		function clear_field(hndl, fdotpath)
			set_field(hndl, fdotpath, []);
		end
		
		% Removes the fieldname corresponding to fdotpath.
		% The handle must be unlocked.
		function delete_field(hndl, fdotpath)
			hndl.check_locked__();
			dotnames = strsplit(fdotpath, '.');
			hndl.deletef__(dotnames);
		end
		
		% Sets all fields to the empty value ([]) but preserving the nested field names.
		function clear_all_fields(hndl)
			fnames = fieldnames(hndl.sdata);
			
			for i = 1 : numel(fnames)
				s = hndl.sdata.(fnames{i});
				if isstruct(s), s = Structxh(s); end
				
				if isa(s, 'Structxh')
					s.clear_all_fields();
					hndl.sdata.(fnames{i}) = s.sdata;
					
				else
					hndl.sdata.(fnames{i}) = [];
				end
			end
		end
		
	end
	
end
%%
