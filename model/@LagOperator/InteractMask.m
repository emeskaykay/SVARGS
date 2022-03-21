% InteractMask
% Logical matrix indicating location of (non-self) interacting components.
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
function bInteract = InteractMask(Op, varargin)
	k = NumVars(Op);
	opts = parse_nvpairs(varargin, {'reduced', 'logical', false});
	bInteract = (any(Op.Phi, 3) | logical(~opts.reduced & Op.Jmat)) & ~logical(eye(k));
end
%%
