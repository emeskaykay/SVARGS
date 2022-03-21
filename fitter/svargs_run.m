% svargs_run
% SVARGS model fit.
%_____
%
% FitParams:
%	An object of type Svargs_params. See Svargs_params.m
%
% dispflag [ default: 1 ]:
% One of:
%	A fileid obtained from a call to fopen - all output is written to this file.
%	or 1 for std output. 
%	or 0 for no output.
%
% OUTPUT:
% mod:
% struct (or struct array) with fields:
% mod.opsparse:
%	LagOperator object with sparse coefficients.
%
% mod.sigma:
%	Noise covariance matrix.
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
function [mod, res, stats] = svargs_run(data, fitparams, dispflag)
	assert( isa(fitparams, 'Svargs_params') );
	assert( isvector(fitparams.LagsAllowed) );
	assert( is_unique(fitparams.LagsAllowed) );
	
	% Set default arguments.
	if (nargin < 3), dispflag = 1; end

	lags = fitparams.LagsAllowed;
	if ismatrix(data)
		fitter = Svargs_fitter(data, fitparams, dispflag);
		[coeffs, res, sigma, intercept] = fitter.Run();
		
		stats = fitter.stats;
		
		addprop(stats, 'intercept');
		stats.intercept = intercept;
		
	elseif ndims(data) == 3
		[k, len, ntrials] = size(data); %#ok<ASGLU>
		fitparams = fitparams.SetLagMultiplier(ntrials);

		data = interleave_trials(data);
		
		% intercept property added to stats during call to svargs_run() aobve.
		[coeffs, sigma, res, stats] = svargs_run(data, fitparams, dispflag);
		[coeffs, ~, res] = relag(coeffs, lags, ntrials, res);
		
	else
		assert(false);
	end
	
	mod.opsparse = LagOperator(lags, coeffs);
	mod.sigma = sigma;
end

% interleave_trials
% Interleave the data. [trial1(:, i), trial2(:, i), .... trialn(:, i) | i = 1,2,3 ...]
%%
function data = interleave_trials(mtdata)
	[k, len, ntrials] = size(mtdata);
	data = permute(mtdata, [1 3 2]);
	data = reshape(data, [k, ntrials*len]);
end

%%
function [C, Lags, res] = relag(C, Lags, ntrials, res)
	Lags = Lags/ntrials;
	
	if (nargin == 4) && (nargout > 2)
		k = size(C,1);
		reslen = size(res,2);
		res = reshape(res, [k, ntrials, reslen/ntrials]);
		res = permute(res, [1 3 2]);
	end
end

%%
