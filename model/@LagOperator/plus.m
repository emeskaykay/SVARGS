% plus | Overloaded addition operator (+).
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
function Opsum = plus(left, right)
	[left, right] = CheckBinaryOp__(left, right, '+');
	
	k = NumVars(left);

	l = left.lags;		r = right.lags;
	ml = length(l);	mr = length(r);
	i = 1;				j = 1;
	S = zeros( k, k, ml+mr );
	Lags = zeros( 1, ml+mr );
	n = 1;
	
	% Note: The lag vectors l(i) and r(j) must be sorted
	% in increasing order for this to work correctly.
	while ( (i <= ml) && (j <= mr) )
		
		if ( l(i) < r(j)  )
			S(:,:,n) = left.Phi(:,:,i);
			Lags(n) = l(i);
			i = i+1;
		elseif ( l(i) > r(j) )
			S(:,:,n) = right.Phi(:,:,j);
			Lags(n) = r(j);
			j = j+1;
		else % l(i) == r(j)
			S(:,:,n) = left.Phi(:,:,i) + right.Phi(:,:,j);
			Lags(n) = l(i);
			i = i+1; j = j+1;
		end
		
		n = n + 1;
	end
	
	while (i <= ml) % r(j) was exhausted.
		S(:,:,n) = left.Phi(:,:,i);
		Lags(n) = l(i);
		i = i+1;
		n = n+1;
	end
	
	while (j <= mr) % l(i) was exhausted.
		S(:,:,n) = right.Phi(:,:,j);
		Lags(n) = r(j);
		j = j+1;
		n = n+1;
	end
	
	% Wipe left over slots.
	S(:,:,n:end) = [];
	Lags(n:end) = [];
	
	I = eye(k);
	Opsum = LagOperator(Lags, S, 'J', I + (left.Jmat-I) + (right.Jmat-I));
	
end
%%
