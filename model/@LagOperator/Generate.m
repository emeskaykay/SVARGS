% Generate
% Usage: Generate(N, y0, wn)
%
% N =		length of data to generate (including initial values).
%
% y0:		(k x p x ntrials) matrix (p initial vectors).
%
% wn:		(k x L x ntrials) array of noise, where L >= (N-p).
%
% Returns:
%	y: (k x N x ntrials) array of VAR data, where k is the number of variables.
%
% NOTE:
% For faster generation, split the length to be generated into multiple
% trials if it is appropriate.
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
function y = Generate(Op, N, y0, wn)
	k = NumVars(Op);
	p = Order(Op);
	ntrials = size(y0, 3);
	assert(size(wn, 3) == ntrials);
	
	assert( all(size(y0,1,2,3) == [k,p,ntrials]) );
	assert( ndims(wn) == ndims(y0) );
	assert( size(wn,1) == k );
	assert( size(wn,2) >= (N-p) );
	
	% Initialize data.
	y = zeros(k, N, ntrials);
	y(:, 1:p, :) = y0;
	
	C = Op.Coeffs(Op.Lags(), 'reduced', true);
	
	% Strech out C as a matrix to do a dot product.
	Cmat = C(:,:);
	
	% Generate data values
	q = numel(Op.lags);
	for t = p+1 : N
		yy = reshape(y(:, t-Op.lags, :), [k*q, ntrials]);
		wt = reshape(wn(:, t-p, :), [k, ntrials]);
		y(:,t, :) = Cmat*yy + wt;
	end
end
%%
