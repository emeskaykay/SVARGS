% [opnew, sr_old, sr_new] = Stabilize(Op, sr_required)
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
function [opnew, sr_old, sr_new] = Stabilize(Op, sr_new, sr_thresh, varargin)
	enforce_modifier(nargout);
	
	opts = parse_nvpairs(varargin, {'sr_old', 'double', []});
	if isempty(opts.sr_old)
		opts.sr_old = model.Operator().SpectralRadius();
	end
	assert(isscalar(opts.sr_old));
	sr_old = opts.sr_old;
	
	if (nargin < 3) || isempty(sr_thresh), sr_thresh = Inf; end
	
	k = Op.NumVars();
	p = Op.Order();
	if (nargin < 2) || isempty(sr_new), sr_new = max(0.99, 1-1/k/p); end
	
	stabilized = (sr_old <= sr_new);
	if stabilized
		opnew = Op;
		sr_new = sr_old;
		disp_x('LagOperator:Stabilize: Operator is already stable upto the required spectral radius. ');
		
	elseif sr_old <= sr_thresh
		sr_ratio = sr_new/sr_old;
		opnew = Op.Reflate(sr_ratio);
		disp_x('LagOperator:Stabilize: Model was stabilized. ');
		
	else
		warning('LagOperator:Stabilize: Spectral radius is greater than stabilization threshold.');
		opnew = Op;
		sr_new = sr_old;
	end
	
end
%%
