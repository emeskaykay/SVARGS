% VariablesCache
% Helps to cache variables conveniently, particularly inside objects.
%
% Due to pass by value semantics in Matlab, data cannot be directly cached
% inside objects in a manner that is transparent to the user. If data is cached
% inside a (non-handle) object property, the function assigning to the property
% becomes a modifier, so that the caller has to copy the modified object back
% into his workspace which can be highly inconvenient to say the least.
%
% The VariablesCache handle would typically be created in the constructor of some
% object and its handle stored as a Transient property (So that it is not saved
% if or when the object is serialized).
%
% Note that all copies of objects containing a cache handle, use the same cache.
% The function Clear() clears the cache for all copies of the object. If an object
% wants to refresh only its own cache variables, it can call the static function
% VariablesCache.New() for example like:
%	hExistingCache = VariablesCache.New(); % Now object has its own empty cache.
%
% Since VariablesCache is implemented using dynamic field names, the names of
% objects to be stored must be valid (possibly dot separated) field names. For
% example "my_func.argumentid.xyz' is a valid name, To store objects that are
% identified by arbitrary data patterns, use the class Cache (see Cache.m), which
% can use arbitrary data id's, but has a longer lookup time.
%
% One can use the VariablesCache class for objects that are the result of fairly
% short calculations because the overhead for the cache lookup equals overhead for
% struct field look up time + overhead of a couple of function calls. The class
% uses Structxh in its implementation.
% (see Structxh.m).
% 
% See also: Cache.m.
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
classdef VariablesCache < handle

	properties (Constant)
		ERR_NAMES = 'VariablesCache: Inputs must all be valid (possibly dot-separated) variable names. '
	end
	
	properties
		hdata;
	end
	
	methods(Static)
		% Return handle to a new Cache. The old cache will be presumably be destroyed
		% by the Java garbage collector if there are no existing referenes  to it.
		function hObj = New(varargin)
			hObj = VariablesCache(varargin{:});
		end
	end
	
	methods
		% varargin: cellarray of property names.
		function hndl = VariablesCache(varargin)
			assert(isempty(varargin) || iscellstr(varargin), VariablesCache.ERR_NAMES);
			numfields = numel(varargin);
			vals = cell(numfields, 1);
			hndl.hdata = Structxh(cell2struct(vals, varargin(:), 1));
		end
		
		% Returns true if there are no values (or only empty values) attached to
		% the existing property names.
		function b = IsEmpty(hndl)
			b = (hndl.hdata.value_count(Inf, true) == 0);
		end
		
		% Get value of property with specified name.
		function val = Get(hndl, fname)
			val = hndl.hdata.get_field(fname);
		end
		
		% Set the value of a property with a specified name. If the property name
		% is not present in the cache, a new property is created and its value set.
		% It would be good programming practice to specify all the property names
		% that might be used in the constructor.
		function Set(hndl, fname, val)
			hndl.hdata.set_field(fname, val);
		end
		
		% True if the property name exists.
		function b = Has(hndl, fname)
			b = hndl.hdata.has_field(fname);
		end
		
		% Get deep copy.
		function cacheobj = GetCopy(This)
			cacheobj = VariablesCache.New();
			cacheobj.hdata = This.hdata.get_copy();
		end
		
		% Clears all property names and values.
		function Clear(hndl)
			hndl.hdata.clear_all_fields();
		end
	end
	
end
%%
