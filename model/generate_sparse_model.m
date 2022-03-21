% [mod, found, sr] = generate_sparse_model(k, p, numcoeffs, csizemin, factorords, sigma, varargin)
% Generates random sparse stable VAR model.
%
% k:
%	Number of variables
%
% p:
%	Scalar value specifying the order of the model.
%
% numcoeffs:
%	Number of coefficients.
%
% csizemin [default: [0, Inf]]:
%	Pair of values, [csizemin, csizemax], specifying minimum and maximum
%	absolute coefficient size. The lower limit cannot be too large: for example
%	a lower limit greater than 0.5 is almost certain to take an impracticable
%	amount of time to find a stable model.
%
% factorords [default = p]:
%	A vector of positive integers. The VAR operator will be found as a product
%	of random VAR operators having the orders specified in factorords. The sum
%	of the integers in factorords must be equal to the order p.
%	When factorords == p (default), the average coefficients size (as a function
%	of lag) has an exponential profile. When factorords is a vector of two or more
%	values, the average coefficient size (as a function of lag) has more complex
%	profile.
%
% varargin:
% Name value pairs
% Name: 'srmax', Value: <double> between 0 and 1. [default: See denude_coeffs()):
%	Upper bound on the admissible spectral radius.
%
% Name: 'sigma', Value: <double> matrix [default: Random noise covariance matrix]:
%	Noise covariance. In case multiple models are expected (for example when runparallel,
%	is true), this can be an array whose each slice are noise-covariance matrices.
%
% Name: 'maxruntime', Value: <double> scalar [default 120 seconds]:
%	Maximum time to run the search.
%
% Name: 'bfastabort', Value: <logical> [default: false]:
%	See denude_coeffs.m
%
% Name: 'nmaxtries', Value: Positive integer [default: 2]
%	Number of times to retry if model is not found within the specified runtime.
%
% Name: 'runparallel', Value: <logical> [default: true]:
%	Runs 'nmaxtries' sparsifications in parallel. In this case upto 'nmaxtries'
%	stable models can be obtained in one run.
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
function [mod, found, sr] = generate_sparse_model(k, p, numcoeffs, csizelims, factorords, varargin)
	if (nargin < 5) || isempty(factorords), factorords = p; end
	if (nargin < 4) || isempty(csizelims), csizelims = [0, Inf]; end
	
	def_maxruntime = max(40, k*p/10 + 2*numcoeffs*csizelims(1));

	opts = parse_nvpairs(varargin, ...
		{'srmax', 'double', []}, ...
		{'sigma', 'double', []}, ...
		{'maxruntime', 'double', def_maxruntime}, ...
		{'bfastabort', 'logical', false}, ...
		{'nmaxtries', 'double', 2}, ...
		{'runparallel', 'logical', false}, ...
		{'verbose', 'logical', false});
	
	srmax = opts.srmax;
	sigma = opts.sigma;
	maxruntime = opts.maxruntime;
	bfastabort = opts.bfastabort;
	nmaxtries = opts.nmaxtries;
	runparallel = opts.runparallel;
	verbose = opts.verbose;
	
	assert(isnumeric(factorords) && isvector(factorords) && all(factorords > 0), ...
		'Factor orders must be vector of positive integers. ');
	assert(sum(factorords) == p, 'Factor orders must sum to order of model. ');
	
	found = false;
	if runparallel
		% Implicit Initialization
		opsparse(nmaxtries,1) = LagOperator();
		found(nmaxtries,1) = false;
		sr(nmaxtries,1) = Inf;
	end
	
	elapsed = 0;
	numfound = 0;
	
	while (numfound == 0) && (elapsed < maxruntime)
		% Note: Strictly, the product of two stable operators (as defined in the product
		% functions below) should be stable, but due to numerical error in multiplication,
		% the product may occassionally fail to be stable.
		opdense = ProductOperatorSeq__(factorords, k, p);
		
		disp_x('Denuding dense operator coefficients. ');
		if runparallel
			poolobj = gcp('nocreate');
			if isempty(poolobj)
				parpool; % create parallel pool
			end
			
			tstart = tic;
			parfor i = 1 : nmaxtries
				[opsparse(i), found(i), sr(i)] = denude_coeffs(opdense.Coeffs(1:p), numcoeffs, csizelims, srmax, bfastabort, verbose);
				if found(i)
					disp(['Found: coeff count: ', num2str(opsparse(i).CoeffCount()), '; Order: ', num2str(opsparse(i).Order())]);
				else
					disp('Not Found');
				end
			end
			
		else
			ntries = 0;
			tstart = tic;
			while ~found && (ntries < nmaxtries)
				[opsparse, found, sr] = denude_coeffs(opdense.Coeffs(1:p), numcoeffs, csizelims, srmax, bfastabort, verbose);
				if found
					disp(['Found: coeff count: ', num2str(opsparse.CoeffCount()), '; Order: ', num2str(opsparse.Order())]);
				else
					disp('Not Found');
				end
				ntries = ntries + 1;
				if ~found && (ntries < nmaxtries)
					disp('Retrying ...');
				end
			end
		end

		numfound = nnz(found);
		elapsed = elapsed + toc(tstart);
	end
	
	fprintf('%u %s\n', numfound, ' stable model(s) found. ');
	
	if numfound == 0
		sr = [];
	else
		opsparse = opsparse(found);
		sr = sr(found);
	end
	
	for i = 1 : numfound
		mod(i).opsparse = opsparse(i); %#ok<AGROW>
		if (size(sigma,1) == 0) || (i > size(sigma,3)) % Matlab reports size([],3) as 1!
			sigma = rand_cov(k,1);
		end
		C0 = opsparse(i).Coeffs(opsparse(i).Lags());
		c = max(abs(C0(:,:)), [], 2);
		c = mean(c);
		mod(i).sigma = c*sigma; %#ok<AGROW>
	end
end

%%
function opx = ProductOperatorSeq__(factorords, k, p)
	progress = Progress_display('Generating VAR operator sequence');
	
	q = factorords(1);
	eiglims = eig_bounds__(k*q, q);
	opx =  var_rand_operator(k, q, eiglims);
	progress.update(1);
	
	for i = 2 : numel(factorords)
		q = factorords(i);
		eiglims = eig_bounds__(k*q, q);
		op = var_rand_operator(k, q, eiglims);
		opx = op + opx - op*opx;
		progress.update(i);
	end
	progress.complete(p);
end

%%
function opx = ProductOperatorTree__(factorsize, k, p)
	progress = Progress_display('Generating VAR operator tree');
	
	p1 = floor(p/2);
	p2 = ceil(p/2);
	lev = 0;
	if p > factorsize
		op1 = ProductOperatorTree__(factorsize, k, p1);
		op2 = ProductOperatorTree__(factorsize, k, p2);
		opx = op1 + op2 - op1*op2;
		
		lev = lev + 1;
		progress.update(2^lev);
		
	else
		eiglims = eig_bounds__(k*p, p);
		opx = var_rand_operator(k, p, eiglims);
		progress.update(2^lev);
	end
	progress.complete(2^lev);
end

%%
function eiglims = eig_bounds__(m, px)
	% Experiment shows that the about 10% to 50% of the eigenvalue magnitudes
	% in a sparse model are close to 0.
	%
	% In addition, as csizemin is increased, the non-zero eigenvalues are pushed
	% towards magnitude 1. We emulate that here to give VARDenudeModel() a better
	% change of finding a sparse model satisfying the minimum size condition.
	%
	% Pushing the magnitudes towards one also leads to a better coefficient spread
	% in the sparsified model - otherwise the coefficients tend to clump together
	% near the minimum specified magnitude.
	na = 0.5*(m - rem(m,2)); % na = Number of complex conjugate eigenvalue pairs.
	eigs_min = rand(na,1);
	eigs_max = eigs_min + rand(na,1).*(1 - eigs_min);
	
	eiglims = [eigs_min, eigs_max].^(1/px);
	eiglims = AdjustEigs__(eiglims);
	if rem(m,2) > 0
		eiglims = [eiglims; [0, 0.01]];
	end
	
	assert(all(abs(eiglims(:)) < 1)); % paranoid check.
end

%%
function Eigs = AdjustEigs__(Eigs)
	r = abs(Eigs);
	rsqr = r.^2;
	rsqradj = icdf('beta', rsqr, 3, 3);
	Eigs = sqrt(rsqradj./rsqr).*Eigs;
end
%%
