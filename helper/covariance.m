% COVARIANCE
% X:
%	(k x N) array of numbers where k is the number of variables,
%	and N is the length of the data.
%
% Y [Optional]:
%	(k x N) array of numbers where k is the number of variables,
%	and N is the length of the data. If Y is not specified, then Y = X.
%
% mx [Optional]:
%	Mean of the x-components.
%
% my [Optional]:
%	Mean of the y-components.
%
% varargin (Optional):
% List of name value pairs:
% 'divisor', d [default: N]: The divisor to use when taking the mean of the squares.
% 'bootstrap', true/false: Use bootstrap to improve estimate.
% 'nsamples', ns [default: 100]: The number of times the bootstrsp is computed.
%
% G:
%	The return value G is the (k x k) covariance Matrix.
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
function G = covariance(X, Y, mx, my, varargin)

	[~, N] = size(X);
	opts = parse_nvpairs(varargin, ...
		{'divisor', 'double', N}, ...
		{'bootstrap', 'logical', false}, ...
		{'nsamples', 'double', 100});
	
	if (nargin < 4), my = []; end
	if (nargin < 3), mx = []; end
	if (nargin < 2), Y = []; end
	
	if ~opts.bootstrap
		if isempty(Y)
			G = covx__(X, mx, opts.divisor);
		else
			assert( all(size(Y) == size(X)) );
			G = covxy__(X, Y, mx, my, opts.divisor);
		end
		
	else
		if isempty(Y)
			G = covx_bs__(X, mx, opts.nsamples, opts.divisor);
		else
			assert( all(size(Y) == size(X)) );
			G = covxy_bs__(X, Y, mx, my, opts.nsamples, opts.divisor);
		end
	end

end

%%
function G = covxy__(X, Y, mx, my, divisor)
	if isempty(my), my = mean(Y,2); end
	if isempty(mx), mx = mean(X,2); end
	G = (X*Y')/divisor - mx*my';
end

function G = covxy_bs__(X, Y, mx, my, nsamples, divisor)
	[k, N] = size(X);
	
	G(k,k) = 0;
	for i = 1 : nsamples
		[Xs, idx] = datasample(X, N, 2);
		Ys = Y(:, idx);
		G = G + covxy__(Xs, Ys, mx, my, divisor);
	end
	G = G/nsamples;
end

%%
function G = covx__(X, mx, divisor)
	if isempty(mx), mx = mean(X, 2); end
	
	% Matlab ensures that (x*x') is hermitian.
	G = (X*X')/divisor - mx*mx';
end

function G = covx_bs__(X, mx, nsamples, divisor)
	[k, N] = size(X);
	
	G(k,k) = 0;
	for i = 1 : nsamples
		Xs = datasample(X, N, 2);
		G = G + covx__(Xs, mx, divisor);
	end
	G = G/nsamples;
end
%%

