% Op = StoreAsCompressed(Op, bstate)
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
function Op = StoreAsCompressed(Op, bstate)
	enforce_modifier(nargout);
	
	Op.b_store_as_compressed = bstate;
	if bstate
		Op.Jmat__ = sparse(Op.Jmat__(:,:));
		Op.Phi___ = sparse(Op.Phi___(:,:));
		
	else
		Op.Jmat__ = full(Op.Jmat__);
		Op.Phi___ = full(Op.Phi___);
		k = size(Op.Phi___, 1);
		Op.Phi___ = reshape(Op.Phi___, k,k,[]);
	end
	
end
%%
