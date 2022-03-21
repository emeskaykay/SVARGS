% H = TransferFunction(Op, Freq)
% Usage: Op.TransferFunction(Freq)
%
% Freq:
%	The set of frequencies at which the transfer function is to be computed.
%	These must be in the half open range [-0.5 : 0.5).
%
% Returns:
% H:
%	The (matrix) transfer function of the linear operator
%		L(y)(t) = e(t)
%	when e(t) is regarded as the input and y(t) the output of the system.
%
%	Here:
%		L(y)(t) = J.y(t) - A1.y(t-1) - A2.y(t-2) - ... - Ap.y(t-p)
%				  = (J - Phi(L))(y(t))
%
%	So if Sigma(w) is the spectrum of e(t), the the spectrum of y(t) is given by:
%		F(w) = H(w)*Sigma(w)*H(w)'
%	where
%		H = inv(T)
%	and
%		T(w) = J - A1.exp(-iw) - A2.exp(-2iw) - ... - Ap.exp(-piw)
%
%	If e(t) is white noise, then Sigma is constant and coincides with the
%	covariance of e(t).
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
function H = TransferFunction(Op, Freq)
	k = NumVars(Op);
	
	m = length(Freq);
	assert(numel(Freq) == m);
	assert( (min(Freq) >= -0.5) && (max(Freq) <= 0.5) );
	Freq = transpose(Freq(:));
	
	I = eye(k);
	nzlags = Op.lags(:);

	if isempty(nzlags)
		H = I(:,:,ones(1,m));
		
	else
		q = length(nzlags);
		
		Jvec = Op.Jmat(:);
		Jrep = Jvec(:, ones(1,m));
		Phi_mat = reshape(Op.Phi, k*k, q);
		
		E = exp(complex(0, -nzlags*2*pi*Freq));
		Hinv_mat = Jrep - Phi_mat*E;
		Hinv = reshape(Hinv_mat, [k,k,m]);
		
		H(k, k, m) = 0; % Implicit Initialization
		
		% *Expensive*. Time complexity is dominated by this loop (90+%).
		for j = 1 : m
			H(:, :, j) = Hinv(:, :, j)\I;
		end
	end
	
end
%%
