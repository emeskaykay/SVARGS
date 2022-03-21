% [ok relerr] = check_relative_error(A, B, reltol, varargin)
%
% A, B:
%	Matrices or vectors to compare. A and B must have the same size.
%
% reltol:
%	Relative error tolerance.
% 
% Name-Value pairs
% Name:
%	errtype,
% Value [default: 'error_in_norm']:
%	String specifying how to reference the relative error. Should be one of
%		'error_in_norm' (=> relerr = ||B-A||/||Ref||)
%		'norm_of_errors'. (=> relerr = ||(B-A)./Ref||
%	where Ref = max(abs(A), abs(B)) + ((A == 0) | (B == 0)); 
%
% Name:
%	fnorm
% Value [default: 2-norm]:
%	One of the following: Either a positive integer OR one of the following
%	strings: 'max', 'euclidean', 'frobenius', 'spectral'. Specifying 'spectral'
%	or 'euclidean' is equivalent to specifying the integer 2.
%
% Name:
%	normtype
% Value [default: 'vector']:
%	One of: 'vector', 'matrix'. If the type is 'matrix', then A and B cannot
%	have dimension greater than 2.
%
% Name:
%	ztol
% Value [default: 0]:
%	Non-negative real number. Values in A or B whose absolute value is less
%	than ztol are set to 0.
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
function [ok, relerr] = check_relative_error(A, B, reltol, varargin)
	assert(isnumeric(A) && isnumeric(B));
	assert(isnumeric(reltol) && isscalar(reltol) && (reltol >= 0));
	assert(all(size(A) == size(B)));
	
	opts = parse_nvpairs(varargin, ...
		{'normtype', 'char', 'vector'}, ...
		{'fnorm', '', 2}, ...
		{'errtype', 'char', 'error_in_norm'}, ...
		{'ztol', 'double', 0});
	if isempty(opts.fnorm), opts.fnorm = 2; end
	assert(isnumeric(opts.fnorm) || ischar(opts.fnorm), ...
		'Norm specification must be a positive integer or character string. ');
	assert(isscalar(opts.ztol) && (opts.ztol >= 0));
	
	normtype = opts.normtype;
	fnorm = opts.fnorm;
	errtype = opts.errtype;
	ztol = opts.ztol;
	
	Ap = abs(A);
	Bp = abs(B);
	if (ztol > 0)
		Ap(Ap < ztol) = 0;
		Bp(Bp < ztol) = 0;
	end
	Ref = max(Ap, Bp);
	Ref = Ref + (Ref == 0);
	
	if strcmpi(errtype, 'norm_of_errors')
		if strcmpi(normtype, 'matrix')
			relerr = matrix_norm__((B-A)./Ref, fnorm);
			
		elseif strcmpi(normtype, 'vector')
			relerr = vector_norm__((B-A)./Ref, fnorm);
			
		else
			assert(false, ['unknown normtype argument: ', normtype]);
		end
		
	elseif strcmpi(errtype, 'error_in_norm')
		if strcmpi(normtype, 'matrix')
			relerr = matrix_norm__(B-A, fnorm)/matrix_norm__(Ref, fnorm);

		elseif strcmpi(normtype, 'vector')
			relerr = vector_norm__(B-A, fnorm)/vector_norm__(Ref, fnorm);
			
		else
			assert(false, ['unknown normtype argument: ', normtype]);
		end
		
	else
		assert(false, ['unknown errtype argument: ', errtype]);
	end
	
	ok = (relerr <= reltol);
end

%%
function a = matrix_norm__(A, fnorm)
	if isnumeric(fnorm)
		p = fnorm;
		a = norm(A,p);
		
	elseif strcmpi(fnorm, 'max')
		a = max(abs(vecc(A)));
		
	elseif strcmpi(fnorm, 'frobenius')
		a = norm(vecc(A), 2);
		
	elseif strcmpi(fnorm, {'euclidean', 'spectral'})
		a = norm(A,2);
		
	else
		assert(false, ['Unknown or unimplemented norm: ' num2str(fnorm)]);
	end
end

%%
function a = vector_norm__(A, fnorm)
	A = vecc(A);
	
	if isnumeric(fnorm)
		p = fnorm;
		a = norm(A,p);
		
	elseif strcmpi(fnorm, 'max')
		a = max(abs(A));
		
	elseif any(strcmpi(fnorm, {'frobenius', 'euclidean', 'spectral'}))
		a = norm(A,2);
		
	else
		assert(false, ['Unknown or unimplemented norm: ' num2str(fnorm)]);
	end
end
%%
