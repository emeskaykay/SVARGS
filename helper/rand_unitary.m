% rand_unitary
% Random unitary n x n matrix. (Random rotation and reflection)
%
% n:
%	length of the required matrix.
%
% varargin:
% Name value pair:
%	name: 'rot', value: true/false [default: false]
%		If 'rot' is true, the result is a random rotation matrix.
%	name: 'uniform_cosines', value: true/false [default: true]
%		If 'uniform_cosines' is true, the direction cosines of the axes are
%		uniformly distributed rather than the angles. If false, it is the
%		angles that are uniformly distributed.
%
% Returns:
% U:
% Unitary or Rotation matrix.
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
function U = rand_unitary(n, varargin)

	opts = parse_nvpairs(varargin, ...
		{'rot', 'logical', false}, ...
		{'uniform_cosines', 'logical', true});
	
	if opts.uniform_cosines
		% Uniformly distributed direction cosines.
		% Reference:
		%	"How to generate random matrices from the classical compact groups",
		%		Francesco Mezzadri, Notices of ACM, Volume 54, Number 5, 2007.
		
		A = randn(n, n); % Normally distributed random matrix.

		[U, R] = qr(A);
		U = U*diag(sign(diag(R)));

	else
		% Uniformly distributed angles.
		% Use a set of random reflections of the identity and a set
		% of random plane rotations to generate the unitary matrix.
		sign__ = 2*(rand(1, n) > 0.5) - 1;
		
		U = diag(sign__);
		
		% I = eye(n);
		T = -pi + 2*pi*rand(n); % Generate random angles beforehand.
		S = sin(T);
		C = cos(T);
		for i = 1 : n
			for j = (i+1) : n
				% R = I;
				% R([i,j], [i,j]) = [c, -s; s, c];
				% U = R*U;
				c = C(i,j);
				s = S(i,j);
				U([i,j], :) = [c, -s; s, c]*U([i,j], :);
			end
		end
	end

	if opts.rot && (det(U) < 0)
		U(:, 1) = -U(:, 1);
	end
	
end
%%
