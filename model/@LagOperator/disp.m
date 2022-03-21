% disp | Overloaded display method.
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
function disp(Op)

	n = numel(Op);
	if n == 1
		for i = 1 : length(Op.lags)
			disp(['Coefficient at lag ', num2str(Op.lags(i)), ' =']);
			disp(Op.Phi(:,:,i));
		end
		
		if Op.has_lag_zero
			disp('Instantaneous multiplier matrix = ');
			disp(Op.Jmat);
		end
		
	else
		for i = 1 : numel(This)
			sk = num2str(NumVars(This(i)));
			so = num2str(Order(This(i)));
			snc = num2str(CoeffCount(This(i)));
			disp(['LagOperator Object <k=', sk, '>, <ord=', so, ...
				'>, <non-zero coeffs=' snc, '>. ']);
		end
	end
	
end
%%
