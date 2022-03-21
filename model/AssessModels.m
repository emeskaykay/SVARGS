% AssessModels
%
% ReferenceModels:
%	Vector of reference structs with fields:
% 	ReferenceModels.opsparse
%	ReferenceModels.sigma
%
% FittedModels:
% Vector of fitted structs with fields:
% 	FittedModels.opsparse
%	FittedModels.sigma
%
% TopCount:
%	If specified will consider only the TopCount Largest coefficients from the fitted
%	models. TopCount can also be a vector of the same length as the first two arguments,
%	in which case, for each i, TopCount(i) will be used to compare ReferenceModels(i)
%	with FittedModels(i).
%	If TopCount == 'refcount', then TopCount(m) is taken to be the number of non-zero
%	coefficients in ReferenceModels(m).
%
% Summary:
% Struct with the following fields:
%	ModelCount: Number of models tested.
%	Count: Total number of coefficients in actual models.
%	Predicted:	Number of locations from actual models present in fitted models.
%	PredictedTop:	Number of 'Count' largest coefficients from actual models present in fitted models.
%	Unpredicted: Number of locations in actual models _not_ present in fitted models
%	Spurious: Number of locations in fitted models not present in actual models.
%	RmsRef: RMS of actual coefficients.
%	RmsUnpredicted: RMS of unpredicted coefficients.
%	RmsError: RMS error of fitted coefficients within actual coefficients.
%	RmsNormalizedError: RMS of relative error of fitted coefficients within actual coefficients.
%	SpuriousWeight: Sqrt of (spurious sum of squares / reference sum of squares).
%	MeanBias: Mean bias of fitted coefficients within actual coefficients.
%	ErrorStdDev: Std Deviation of RMSError.
%	OrderError: Normalized sum, over spurious lags, of the lagwise maximal coefficient.
%	Correlation: Average correlation between coefficients of fitted and actual models.
%	NoiseDiff: Average norm of relative differences in noise terms between fitted models actual models.
%	
% FieldText:
%	Printable Text for each field name. Same field names as Summary.
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
function [Summary, FieldText] = AssessModels(ReferenceModels, FittedModels, TopCount)
	
	nmods = length(ReferenceModels);
	
	if nargin < 3, TopCount = []; end
	if ischar(TopCount)
		assert(strcmpi(TopCount, 'refcount'));
		TopCount = zeros(1, nmods);
		for m = 1 : nmods
			TopCount(m) = ReferenceModels(m).opsparse.CoeffCount();
		end
	end
	assert(isempty(TopCount) || isvector(TopCount));
	assert(isempty(TopCount) || isnumeric(TopCount));
	if isscalar(TopCount), TopCount = TopCount(:, ones(1,nmods)); end
	
	correl(nmods, 1) = 0;
	ref_detnoise(nmods, 1) = 0;
	fit_detnoise(nmods, 1) = 0;
	noise_diffs(nmods, 1) = 0;
	
	refvec = [];
	fitvec = [];
	fitvec_top = [];
	maxcoeffs_at_spurious_lags = [0;0];
	for m = 1 : nmods
		refmod = ReferenceModels(m);
		fitmod = FittedModels(m);
		if isempty(fitmod.opsparse), continue; end
		
		fit_detnoise(m) = 0;
		[~, p] = chol(fitmod.sigma);
		if p == 0
			fit_detnoise(m) = log_determinant(fitmod.sigma);
		end
		ref_detnoise(m) = log_determinant(refmod.sigma);
		
		[~, noise_diffs(m)] = check_relative_error(refmod.sigma, fitmod.sigma, Inf);
		
		ordref = refmod.opsparse.Order();
		ordfit = fitmod.opsparse.Order();
		ord = max(ordref, ordfit);
		refcoeffs = refmod.opsparse.Coeffs(1:ord);
		fitcoeffs = fitmod.opsparse.Coeffs(1:ord);
		if ~isempty(TopCount)
			fitcoeffs_top = KeepLargest(fitcoeffs, TopCount(m), @abs);
		else
			fitcoeffs_top = fitcoeffs;
		end
		
		lagsref = refmod.opsparse.Lags();
		lagsfit = fitmod.opsparse.Lags();
		spurious_lags = setdiff(lagsfit, lagsref);
		fitcoeffs_spurlags = max(max(abs(fitcoeffs(:,:,spurious_lags)), [], 1), [], 2);
		maxcoeffs_at_spurious_lags = [maxcoeffs_at_spurious_lags; vecc(fitcoeffs_spurlags)]; %#ok<AGROW>
		
		bnz = ( (refcoeffs(:) ~= 0) | (fitcoeffs(:) ~= 0) );
		refcoeffs = vecc(refcoeffs(bnz));
		fitcoeffs = vecc(fitcoeffs(bnz));
		fitcoeffs_top = vecc(fitcoeffs_top(bnz));
		
		% Need to use vecc, otherwise fails for k = 1 inconsistent Matlab behaviour.
		correl(m) = corr(vecc(refcoeffs(refcoeffs ~= 0)), vecc(fitcoeffs(refcoeffs ~= 0)));
		
		refvec = [refvec; refcoeffs]; %#ok<AGROW>
		fitvec = [fitvec; fitcoeffs]; %#ok<AGROW>
		fitvec_top = [fitvec_top; fitcoeffs_top]; %#ok<AGROW>
	end

	diffvec_top = (fitvec_top - refvec);
	
	bref = (refvec ~= 0);
	bpred = (refvec & fitvec_top);
	bun = (refvec & ~fitvec_top);
	bspur = (~refvec & fitvec);
	
	diffvec_top_normalized = diffvec_top(bref)./refvec(bref);
	fitvec_spur = fitvec(bspur);

	Summary.ModelCount = nmods;
	Summary.Count = nnz(bref);
	Summary.Predicted = nnz(bpred);
	Summary.Unpredicted = nnz(bun);
	Summary.Spurious = nnz(bspur);
	Summary.SpuriousWeight = sqrt( sum(fitvec_spur.*fitvec_spur)/sum(refvec(bref).*refvec(bref)) );
	Summary.RmsRef = sqrt( mean( refvec(bref).*refvec(bref) ) );
	Summary.RmsUnpredicted = sqrt( mean( refvec(bun).*refvec(bun) ) );
	Summary.RmsError = sqrt( mean( diffvec_top(bref).*diffvec_top(bref) ) );
	Summary.RmsNormalizedError = sqrt( mean( diffvec_top_normalized.*diffvec_top_normalized ) );
	Summary.MeanNormalizedBias = mean( diffvec_top_normalized );
	Summary.NormalizedErrorStdDev = sqrt(Summary.RmsNormalizedError^2 - Summary.MeanNormalizedBias^2);
	Summary.OrderError = sum(maxcoeffs_at_spurious_lags)/nmods/sqrt(mean(fitvec(bref).*fitvec(bref)));
	Summary.Correlation = mean(correl);
	Summary.NoiseDiff = mean(noise_diffs);
	
	if nnz(bref) == 0
		Summary.RmsRef = 0;
		Summary.RmsError = 0;
		Summary.MeanNormalizedBias = 0;
		Summary.NormalizedErrorStdDev = 0;
	end
	
	if nnz(bun) == 0
		Summary.Unpredicted = 0;
	end
	
	if nnz(bspur) == 0
		Summary.SpuriousWeight = 0;
	end
	
	FieldText.ModelCount = 'Number of models';
	FieldText.Count = 'Coefficient count in true model';
	FieldText.Predicted = 'Predicted count';
	FieldText.Unpredicted = 'Unpredicted count';
	FieldText.Spurious = 'Spurious count';
	FieldText.RmsRef = 'RMS of coefficient values in true model';
	FieldText.RmsUnpredicted = 'RMS of values of unpredicted coefficients';
	FieldText.RmsError = 'RMS Error';
	FieldText.RmsNormalizedError = 'RMS of normalized Error';
	FieldText.SpuriousWeight = 'Sqrt of (spurious sum of squares / reference sum of squares).';
	FieldText.MeanNormalizedBias = 'Normalized bias in mean over all coefficients. ';
	FieldText.NormalizedErrorStdDev = 'Std. Dev. of normalized error';
	FieldText.OrderError = 'Mean absolute error in order';
	FieldText.Correlation = 'Correlation between predicted and actual coefficients. ';
	FieldText.NoiseDiff = 'Normalized RMS of difference between computed and actual noise covariance matrix. ';
	
end

%%
