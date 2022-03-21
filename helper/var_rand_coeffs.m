% var_rand_coeffs
%
% Generates random (real) coefficients of a k variable VAR(p) process,
% given the eigenvalues of the process. This can be used, for example
% to generate VAR coefficients while ensuring that the process is
% stable (By choosing the eigenvalues inside the unit circle).
% C:				Return value is an array of coefficient matrices.
%
% fwd_eigs:		Eigenvalues of the (forward) characteristic polynomial.
%
% k:				Number of variables in the process.
%
% p:				Order of the process.
%
% The input arguments must satisfy k*p = (number of eigenvalues).
%
% The function works by randomly generating a certain part of each
% eigenvector associated with the corresponding AR(1) process. The
% remaining part of each eigenvector is determined by this part.
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
function C = var_rand_coeffs(k, p, fwd_eigs)
	assert( (k > 0) && (p > 0) );
	assert( length(fwd_eigs) == k*p, 'Number of eigenvalues must be (nVars * Order). ' );
	assert( numel(fwd_eigs) == k*p, 'Eigenvalue argument must be vector. ' );			
	assert( is_unique(fwd_eigs), 'Eigenvalues must be distinct. ');

	fwd_eigs = vecr(fwd_eigs); % Make into row vector in case it is column.
	eigs_cx1 = sort( fwd_eigs(imag(fwd_eigs) > 0) );
	eigs_cx2 = sort( fwd_eigs(imag(fwd_eigs) < 0) );
	eigs_r = fwd_eigs(imag(fwd_eigs) == 0);
	assert( all(eigs_cx2 == conj(eigs_cx1)), 'Complex eigenvalues must be conjugate pairs. ');

	n = length(eigs_r);
	m = length(eigs_cx1);

	blockdiag = zeros(k*p, k*p);
	eigmat = zeros(k*p, k*p);

	if (m > 0)
		rs = zeros(k*p, 2);
		rsp1 = rand_sphere(k, ones(1,m), ones(1,m));
		rsp2 = rand_sphere(k, ones(1,m), ones(1,m));
		for c = 1:m
			a = real(eigs_cx1(c));
			b = imag(eigs_cx1(c));
			for j = 1 : p
				rs((j-1)*k+1 : j*k, :) = [rsp1(:,c) rsp2(:,c)]*[a b; -b a]^(p-j);
			end
			h = 2*c-1;
			blockdiag(h:h+1, h:h+1) = [a b; -b a];
			eigmat(:, h:h+1) = rs;
		end
	end

	if (n > 0)
		Q = zeros(k*p,n);
		ll = repmat(eigs_r, [k 1]);
		qq = rand_sphere(k, ones(1,n), ones(1,n));
		for j = 1:p
			Q((j-1)*k+1 : j*k, :) = qq.*(ll.^(p-j));
		end
		blockdiag(2*m+1:end, 2*m+1:end) = diag(eigs_r);
		eigmat(:, 2*m+1:end) = Q;
	end

	A = eigmat(1:k, :)*blockdiag/eigmat; % k x kp matrix.
	C = reshape(A, [k,k,p]); % Make the coefficients from A.

end
%%
