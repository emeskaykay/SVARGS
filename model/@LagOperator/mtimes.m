% mtimes | Overloaded multiplication (*) operator.
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
function OpProd = mtimes(left, right)
	[left, right] = CheckBinaryOp__(left, right, '*');
	
	if isnumeric(left)
		Phinew = reshape(left*right.Phi(:,:), size(right.Phi));
		
		J = Op.Jmat;
		I = diag(diag(J));
		OpProd = LagOperator(right.lags, Phinew, 'J', I + left*(J-I));
		
	else
		OpProd = Compose__(left, right);
	end
end
%%
