% Coeffs
% Usage: Coeffs(lags)
%
% lags:
%	Vector of positive lags containing the lags at which the coefficient
%	are required (This must be specified).
%
% varargin:
% Name-Value Pairs;
%	Name: 'reduced', Value: true/false. [default: false]
%	If reduced is true, the returned array C are the coefficents multiplied by
%	the inverse of J (the coefficient of Y(t) and the returned value of J is the
%	identity.
%
% Return values:
% C:
%	A (kxkxn) array, where n is the length of the 'lags' vector. Each
%	slice of the array is a matrix coefficient. Only the coefficients
%	corresponding to the lags specified in the input argument (lags)
%	are returned. If a particular specified lag is not present, the
%	corresponding slice in the return array is the zero matrix. 
%
% nzix:
%	Vector containing those indices for which the coefficients are
%	non-zero. Thus C(nzix) gives the non-zero coefficients of C.
%
% J:
%	The coefficient of Y(t). [ J.Y(t) = C1.Y(t-1) + C2.Y(t-1) ... ].
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
function [C, nzix, J] = Coeffs(Op, lags, varargin)
	k = NumVars(Op);
	if isempty(lags)
		C = zeros(k,k,0);
		nzix = [];
		J = Op.Jmat;
		return;
	end
	
	assert(isvector(lags) && is_integral(lags, [1, Inf]));
	
	opts = parse_nvpairs(varargin, {'reduced', 'logical', false});
	
	if isequal(lags, Op.lags)
		C = Op.Phi;
		nzix = 1 : numel(lags);
		
	else
		C(k, k, length(lags)) = 0; % Implicit Initialization

		% The index vectors returned by intersect() satisfy:-
		%	Op.lags(ixop) == lags(nzix) == (common lags).
		[~, ixop, nzix] = intersect(Op.lags, lags);
		
		% Unbelievably expensive operation in matlab.
		C(:,:,nzix) = Op.Phi(:,:,ixop);
	end
	
	J = Op.Jmat;
	if opts.reduced
		C = J\C(:,:);
		C = reshape(C, k,k,[]);
		if nargout > 2
			J = eye(k);
		end
	end
end
%%
