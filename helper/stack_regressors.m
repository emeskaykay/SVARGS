% stack_regressors
% Lag specified columns of the transposed time series data and stack them to
% make a matrix.
%
% data:
%	(k x N) numeric matrix where k is the number of variables and N is the number of time points.
%
% hregcache:
%	handle to a VariableCache object.
%
% btimes:
%	Numeric vector of time indices. After regressors are stacked, only rows
%	corresponding to ixtimes are returned. If this is a scalar (btimes = istart),
%	then rows (istart:end) are returned.
%
% lags:
%	Vector of lags. Columns of Y (the transposed time series) are shifted by the
%	lags specified in this vector. Elements of lags must be integers, positive,
%	negative or zero.
%
% blocs:
%	(k x q) Boolean mask, where q is the number of lags specified in 'lags'. The
%	i-th column of blocs indicates the variable indices that are to be used to form
%	the regressors that are lagged by lag(i).
%
%	If blocs == 'all', then all variables are stacked.
% 
%	If blocs == 'cache_all', then all variables are cached provided that all of the
%	required lags are found in the cache. If one or more of the required lags are not
%	found in the cache, a new stack of regressors is created using the lags specified
%	in the 'lags' variable and the result is cached, overwriting any previously existing
%	regressors in the cache.
%
% baddones:
%	If true, a column of ones is added to the last column of the set of regressors X.
%
% Returns:
% X:
%	The set of regressors formed by lagging columns of transpose(Y).
% 
% Y:
%	Transpose of the time series data.
%
% NOTE:
%	The columns in X are in increasing order of lag.
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
function [X, Y] = stack_regressors(data, hregcache, btimes, lags, blocs, baddones)

	if 0 == exist('baddones', 'var'), baddones = false; end
	assert(islogical(baddones) && isscalar(baddones));
	
	len = size(data, 2);
	if isscalar(btimes)
		assert(isnumeric(btimes));
		istart = btimes;
		assert((istart >= 1) && (istart < len));
		btimes = [false(1, istart-1), true(1, len-istart+1)];
	else
		assert(isvector(btimes) && islogical(btimes));
		assert(numel(btimes) == len);
	end

	% If regressors at the required lags are cached, get them.
	X = get_cached_regressors(hregcache, lags, blocs);
	
	if ~isempty(X)
		Y = transpose(data);

	elseif ischar(blocs)
		if strcmpi(blocs, 'all')
			[X, Y] = stack_all(data, hregcache, lags, false);
		elseif strcmpi(blocs, 'cache_all')
			[X, Y] = stack_all(data, hregcache, lags, true);
		else
			assert(false, 'Invalid ''blocs'' argument. ');
		end
		
	else
		% Remove lags with no coefficients.
		bnzlags = any(blocs, 1);
		blocs = blocs(:, bnzlags);
		lags = lags(bnzlags);
		
		% Stack chosen regressors.
		[X, Y] = stack_choice(data, lags, blocs);
	end
	
	X = X(btimes, :);
	Y = Y(btimes, :);
	
	if baddones
		X(:, end+1) = ones(nnz(btimes), 1);
	end
end

%%
function X = get_cached_regressors(hregcache, lags, blocs)

	Xall = hregcache.Get('Regressors');
	
	if isempty(Xall)
		X = [];
		
	else
		lags_Xall = hregcache.Get('RegressorLags');
		[~, ~, ixlags] = SetIntersection(lags, lags_Xall);
		if (numel(ixlags) ~= numel(lags))
			X = [];

		elseif isequal(blocs, 'cache_all')
			X = Xall(:,:,ixlags);
			n = size(X,1);
			X = reshape(X, n, []);

		else
			X = Xall(:,:,ixlags);
			X = X(:, blocs(:));
			
		end
	end
	
end

%%
function [X, Y] = stack_choice(data, lags, blocs)
	mperlag = sum(blocs, 1);
	mtotal = sum(mperlag);
	[~, n] = size(data);
	
	Y = transpose(data);  % Response vector.
	
	X = zeros(n, mtotal);
	if mtotal > 0
		% Stack up regressors lag by lag.
		last = 0;
		q = length(lags);
		for i = 1 : q
			bchoice = blocs(:, i);
			h = lags(i);
			if h >= 0
				X(h+1:end, last+1:last+mperlag(i)) = Y(1:end-h, bchoice);
				
			elseif h < 0
				g = -h;
				X(1:end-g, last+1:last+mperlag(i)) = Y(g+1:end, bchoice);				
			end
			last = last + mperlag(i);
		end
	end
	
end

function [Xall, Y] = stack_all(data, hregcache, lags, bcache)
	[k, n] = size(data);

	Y = transpose(data); % Response vector.
	
	q = length(lags);
	Xall = zeros(n, k, q);
	for i = 1 : q
		h = lags(i);
		if h >= 0
			Xall(h+1:end, :, i) = Y(1:end-h, :);
			
		elseif h < 0
			g = -h;
			Xall(1:end-g, :, i) = Y(g+1:end, :);
		end
	end
	
	if bcache
		if isempty(Xall)
			disp('Warning: Storing empty regressor array. ');
		end
		hregcache.Set('Regressors', Xall);
		hregcache.Set('RegressorLags', lags);
	end
	
	Xall = reshape(Xall, [n, k*q]);
end

%%
