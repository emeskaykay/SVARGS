% Plot
% Plots a graph of the two specified variables.
%
% xName:
%	Variable to plot on the x-axis. One of:
%		'modelid'
%		'complexity'
%
% yName:
%	Variable to plot on the y-axis. One of:
%		'complexity'
%		'cost'
%		'quality'
%
% xName must not be the same as yName.
%
% hCriterion:
%	If yName is 'quality', then the handle to the criterion method must be specified.
%
% graphObj [optional]:
%	The Graph object on which to plot.
%
% varargin:
%	Matlab parameters for the plot.
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
function [graphObj, xvals, yvals] = Plot(This, xName, yName, hCriterion, graphObj, varargin)

	if ( (nargin < 5) || isempty(graphObj) || ~graphObj.IsValid() )
		graphObj = Graph;
	end

	assert( ~strcmpi(xName, yName) );
	
	switch lower(xName)
		case 'modelid'
			ids = ModelIds(This);
			xvals = ids;
		case 'complexity'
			[ids, xvals] = ModelIds(This, 'sort_complexity');
		otherwise
			error('MOS:Plot:Invalid variable for x-axis. ');
	end
	
	critName = '';
	switch lower(yName)
		case 'complexity'
			yvals = Complexity(This, ids);
		case 'cost'
			yvals = Cost(This, ids);
		case 'quality'
			yvals = Quality(This, hCriterion, ids);
			critName = [ '(', func2str(hCriterion), ' value)' ];
		otherwise
			error('MOS:Plot:Invalid variable for y-axis. ');
	end
	
	bInf = isinf(yvals);
	bMinusInf = bInf & (yvals < 0);
	bPlusInf = bInf & (yvals > 0);
	
	if nnz(bMinusInf) > 0
		if nnz(bMinusInf) == nnz(~bPlusInf)
			yvals(bMinusInf) = -eps;
		else
			yvals(bMinusInf) = min(yvals(~bMinusInf));
		end
	end
	
	if nnz(bPlusInf) > 0
		if nnz(bPlusInf) == nnz(~bMinusInf)
			yvals(bPlusInf) = eps;
		else
			yvals(bPlusInf) = max(yvals(~bPlusInf));
		end
	end
	
	info = PlotInfo(xvals, yvals, [], 'line');
	info.yLabel = [yName, critName];
	info.xLabel = xName;
	
	AxisTitle = info.axisTitle();
	
	graphObj.HoldState(1, 'on');
	bAddStars = false;
	if nnz(bMinusInf) > 0
		bAddStars = true;
		x_minus_inf = xvals(bMinusInf);
		y_minus_inf = yvals(bMinusInf);
		info_minus_inf = PlotInfo(x_minus_inf, y_minus_inf, [], 'line');
		AxisTitle = [AxisTitle, {'Green stars indicate minus infinity (best possible quality). '}];
		graphObj.Plot(info_minus_inf, '*', 'Color', 'green');
	end
	
	if nnz(bPlusInf) > 0
		bAddStars = true;
		x_inf = xvals(bPlusInf);
		y_inf = yvals(bPlusInf);
		info_inf = PlotInfo(x_inf, y_inf, [], 'line');
		AxisTitle = [AxisTitle, {'Red stars indicate infinity (worst possible quality). '}];
		graphObj.Plot(info_inf, '*', 'Color', 'red');
	end
	
	graphObj.Plot(info, varargin{:});
	graphObj.HoldState(1, 'off');
	
	if bAddStars
		AxisTitle = [AxisTitle, {'These values were not in the available normalized floating point range. '}];
		graphObj.AxisTitle(AxisTitle, 'Color', 'blue');
	end

end
%%
