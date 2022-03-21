% Class MOS (Model Order Selector)
%
% N:
%	The length of the data sample that is used to estimate the model (The original
%	sample length minus the length of any presample data or initial values).
%
% k:
%	The dimension (number of variables) of the time series.
%
% M:
%	Maximum number of parameters (regressors). [Only required if sbc_x criterion is used].
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
classdef MOS < handle

	properties (Access = 'private')
		idrecmap = [];
	end
	
	properties (SetAccess = 'private')
		N;
		k;
		M;
	end
	
	methods (Static, Access = 'private')
		irmap = CreateMap();
	end
	
	methods
		function This = MOS(N, k, M)
			if nargin > 0
				This.N = N;
				This.k = k;
				This.M = M;
				This.idrecmap = CreateMap();
			end
		end
		
		N = SampleLength(This);
		c = ModelCount(This);
		
		[ids, vals] = ModelIds(This, sortby, hcrit)
		mvec = Complexity(This, ids);
		svec = Cost(This, ids);
		qvals = Quality(This, hcrit, ids);
		[idbest, valbest] = Best(This, hcrit, ids, mintype, varargin)
		
		[graphObj, xvals, yvals] = Plot(This, xName, yName, hCriterion, graphObj, varargin)
		
		% Modifiers
		Insert(This, ids, complexities, costs, owwarn)
		Remove(This, ids);
		Merge(This, sright);
		
	end
	
end
%%
