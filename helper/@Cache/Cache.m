% Cache
% Helps to cache data conveniently, particularly inside objects.
%
% Due to pass by value semantics in Matlab, data cannot be directly cached
% inside objects in a manner that is transparent to the user. If data is cached
% inside a (non-handle) object property, the function assigning to the property
% becomes a modifier, so that the caller has to copy the modified object back
% into his workspace which can be highly inconvenient to say the least.
%
% The Cache handle would typically be created in the constructor of some object
% and its handle stored as a Transient property (So that it is not saved if or
% when the object is serialized).
%
% Note that all copies of object that containing a cache handle, use the same
% cache. The function Clear() clears the cache for all copies of the object. If
% an object wants to refresh only its own cache variables, it can call the static
% function Cache.New() for example like:
%	hExistingCache = Cache.New(); % Now object has its own empty cache.
%
% Note that Data Id's in the cache need not be valid be variable names. They can
% be any arbitrary strings or arbitrary data patterns that can be compared using
% isequal(). If only variables need to be cached, use the class VariablesCache
% instead, which is more efficient.
%
% It may not pay to use this Cache class for objects that are the results of very
% short calculations because the cache lookup time may outweigh the calculation
% time. For objects that are the results of short calculations, use the class
% VariablesCache instead.
% 
% See VariablesCache.m.
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
%% TODO: Use Containers:map to implement storage?
classdef Cache < handle
	
	properties (Access = 'private')
		DataStorage = {};
		DataIdCells = {};
	end

	methods (Static)
		function hObj = New()
			hObj = Cache();
		end
	end
	
	methods (Access = 'public')
		
		% Get a stored value by its DataID.
		function v = Get(this, DataID)
			v = [];

			k = 1;
			idcount = length(this.DataIdCells);
			while (k <= idcount) && ~isequal(DataID, this.DataIdCells{k})
				k = k + 1;
			end

			if ~isempty(this.DataStorage) && (k <= idcount)
				v = this.DataStorage{k};
			end
		end

		% Associate a data object with an id and store it.
		% The DataID can be any object type that can be compared with isequal().
		function Set(this, DataID, Data)
			idcount = length(this.DataIdCells);
			k = 1;
			while (k <= idcount) && ~isequal(DataID, this.DataIdCells{k})
				k = k + 1;
			end
			if k <= idcount
				this.DataStorage{k} = Data;
			else
				this.DataIdCells{end+1} = DataID;
				this.DataStorage{end+1} = Data;
			end
		end
		
		function b = Has(this, DataID)
			k = 1;
			idcount = length(this.DataIdCells);
			while (k <= idcount) && ~isequal(DataID, this.DataIdCells{k})
				k = k + 1;
			end
			b = (k <= idcount);
		end
		
		% Get deep copy.
		function cacheobj = GetCopy(This)
			cacheobj = Cache.New();
			cacheobj.DataStorage = This.DataStorage;
			cacheobj.DataIdCells = This.DataIdCells;
		end
		
		function Clear(this)
			this.DataStorage = {};
			this.DataIdCells = {};
		end
		
	end

end
%%
