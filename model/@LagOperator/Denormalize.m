% This = Denormalize(This, lambda)
%
% Useful when fitting normalized time series. The resulting lag operator
% then needs to be denormalized to get the lag operator for the original
% data.
%
% Note that this does not work if the operator is the result of fitting
% a trialwise-normalized multitrial series. The multitrial series must
% have been ensemble-normalized.
%
% lambda:
%	The vector of normalizing factors that were used to divide the
%	original time series.
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
function Op = Denormalize(Op, lambda)
	Op.Phi___ = bsxfun(@times, vecc(lambda), Op.Phi___);
	Op.Phi___ = bsxfun(@rdivide, Op.Phi___, vecr(lambda));
	
	% Can probably adjust cache variables, but for
	% the moment get a handle to a fresh cache.
	Op = Op.CacheNew();
end
%%
