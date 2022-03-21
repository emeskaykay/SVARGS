% val = llh(mosobj, logdetsigma, ~)
% Calculates the Log-Likelihood from the log absolute error determinant.
% 
% logdetsigma:
%	log determinant of error covariance. This can be a vector.
%
% L = -0.5*N*( log(det(sigma)) - k*(1 + log(2*pi) )
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
function val = llh(mosobj, logdetsigma, ~)

	k = mosobj.k;
	N = mosobj.N;
	
	val = -0.5*N*( logdetsigma - k*(1 + log(2*pi)) );
	
end
%%
