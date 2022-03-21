% GIVENS
% givens(t)
%	t is the rotation angle.
% givens(a, b)
%	The rotation angle is such that the vector [a, b] is rotated so that the second coordinate is 0.
%		ie: [a, b] ->(t) [cos(t)*a - sin(t)*b, 0];
%
% givens.Block:
%	2x2 rotation matrix in the plane through the angle (t) specified in the constructor.
%
% givens.Rotate(U,i,j):
%	Returns the [i,j] x [i,j] block of the matrix U rotated (clockwise) through the angle (t)
%	specified in the constructor. 
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
classdef givens
	
	properties (SetAccess = private)
		% c == cos(t), s = sin(t), r = sqrt(a*a+b*b)
		c, s, r;
	end
	
	methods
		function This = givens(a, b)
			if (nargin == 1)
				if (a ~= 0)
					This.c = cos(a);
					This.s = sin(a);
				else
					This.c = 1;
					This.s = 0;
				end
				This.r = 1;
			else
				% Effective computation is:
				% This.r = norm([a;b]);
				% This.c = a/This.r;
				% This.s = b/This.r;
				% but we rearrange it below to minimize potential for overflow.
				sgn_a = sign(a);
				sgn_b = sign(b);
				abs_a = a*sgn_a;
				abs_b = b*sgn_b;
				if (abs_a ~= 0) && (abs_b ~= 0)
					if abs_b > abs_a
						t = a/b;
						d = sqrt(1 + t*t);
						This.s = sgn_b/d;
						This.c = t*This.s;
						This.r = abs_b*d;
					else
						t = b/a;
						d = sqrt(1 + t*t);
						This.c = sgn_a/d;
						This.s = t*This.c;
						This.r = abs_a*d;
					end
				elseif abs_b == 0
 					This.r = a;
 					This.c = 1;
 					This.s = 0;
				else
 					This.r = b;
 					This.c = 0;
 					This.s = 1;
				end
			end
		end
		
		% [cos(t), sin(t)] pair.
		function [ct, st] = Pair(This)
			ct = This.c;
			st = This.s;
		end
		
		% Clockwise rotation matrix [c, s; -s, c].
		% * Premultiply by G to rotate column clockwise (givens rotation).
		% * Premultiply by transpose(G) to rotate column anti-clockwise.
		% * The givens (clockwise) rotation zeros the second element in the column vector [a; b].
		function G = Block(This)
			G = [This.c, This.s; ...
					-This.s, This.c];
		end
		
		% Given rotation on U([i,j], [i,j]).
		function u22 = Rotate(This, U, i, j)
			assert(i ~= j);
			if (i > j), [i, j] = Swap(i, j); end
			u22 = Block(This)*U([i,j], [i,j]);
		end
		
	end
end
%%
