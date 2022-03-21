% plot
% Plot and compare coefficients.
%
% mLags [default: all lags upto order]:
%	The lags of the coefficients that are to be plotted.
%
% bInteracts [default: all interactions]:
%	(k x k) logical matrix indicating the (i,j) interacts to be plotted.
%
% g [Default: Creates a new graph object]
%	Graph object where the plot is to take place.(See Graph.m).
%
% gParams.xCoords [Default: Same values as the 'mLags' argument].
%	An increasing sequence of values each giving the x-coordinate
%	of the corresponding lags specified in mLags.
%
% gParams.bPlotMask [Default: all ones]
%	A (k x k) logical matrix indicating which (row, column) pairs
%	of the coefficients to plot.
%
% gParams.rcVec [default: determined by the best fit]
%	Two element vector [nrows, ncols] giving the number of subplots,
%	vertically and horizontally respectively, in the graph.
%
% gParams.bHold [default: false]
%	If the function is called once with gparams.bhold = false, then
%	subsequent calls to _any_ LagOperator object with gparams.bhold
%	set to true, will overlay the coefficients (suitably offset) on
%	the same graph (unless the original graph has been deleted).
%
% gParams.yLabels [default: Coeffs(row # of coeff, col # of coeff)]
%	Cell vector of y-axis labels, in column-first order, for each
%	(row, column) that is to be plotted.
%
% gParams.gCount, gParams.barWidth:
%	Used by the implementation.
%	
% varargin:
%	Args to pass to the 'bar' function
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
function g = plot(Oper, mLags, bInteracts, g, gParams, varargin)

	k = NumVars(Oper);
	if (nargin < 5) || isempty(gParams), gParams = []; end
	if (nargin < 4) || isempty(g), g = []; end
	if (nargin < 3) || isempty(bInteracts), bInteracts = true(k); end
	
	[mLags g gParams varargin] = GetParams__(Oper, mLags, bInteracts, g, gParams, varargin{:});
	
	% If there is only a single lag (or no lags!), add zero coefficients
	% so that MultiAxesPlot knows to plot a 3-dimensional array.
	if length(mLags) == 1
		mLags = [mLags, mLags(end)+1];
	elseif isempty(mLags)
		mLags = [1, 2];
	end
	
	g = MultiAxesPlot(Coeffs(Oper, mLags), mLags, g, gParams, varargin{:});
	g.ScaleAxes('tight');
	g.ScaleAxes([], 'quilt');
	
end

%%
function [mLags g gParams varargin] = GetParams__(Oper, mLags, bInteracts, g, gParams, varargin)
	Ignore(Oper);
	
	gParams.PlotType = 'bar';
	
	if ~isfield(gParams, 'yLabels') || isempty(gParams.yLabels)
		gParams.yLabels = 'Coeffs';
	end
	
	if ~isfield(gParams, 'xLabel') || isempty(gParams.xLabel)
		gParams.xLabel = 'Lag';
	end
	
	if ~isfield(gParams, 'bPlotMask') || isempty(gParams.bPlotMask)
		gParams.bPlotMask = bInteracts;
	else
		gParams.bPlotMask = (gParams.bPlotMask & bInteracts);
	end
end
%%
