% LagOperator
% Object representing a VAR(p) lag operator:
%	J.Y(t) = C1.y(t-1) + C2.y(t-2) + ... + Cp.y(t-p)								[1]
% where J = I - C0 where the diagonal of C0 is zero and I is the identity matrix.
%
% The lag operator can be interpreted in more than one way depending on what
% observations are available at time t. If no observations are available at
% time t, then it is effectively equivalent to
%	Y(t)	= inv(J).[C1.y(t-1) + C2.y(t-2) + ... + Cp.y(t-p)]						[2]
%
% Using the notation L(Y(t)) = Y(t-1), this can be written as:
%	Y(t)	= Phi(L)(Y(t)), where Phi(L) = C1.L + C2.L^2 + ... + Cp.L^p				[3]
%
% On the other hand we can also interpret the equation as:
%	Y(t) = C0.Y(t) + C1.y(t-1) + C2.y(t-2) + ... + Cp.y(t-p) = Phix(L)(Y(t))
% where C0 = I - J.
% Such an interpretation with _all_ the instantaneous terms on the right seems
% unlikely to be practically useful. However, if _some_ of the instantaneous
% values are known, then the equation:
%	(J + C0x).Y(t) = C0x.Y(t) + C1.y(t-1) + C2.y(t-2) + ... + Cp.y(t-p)	[4]
% is relevant. Here C0x is C0 with the columns corresponding to the unknown
% components of Y removed.
%
% For the purpose of propagating Y indefinitely, mostly [2] is useful. Therefore
% methods for this class that require an operator interpretation will use the
% inverted form in [2].
%
% Constructor
% Op = LagOperator(Lags, Coeffs, varargin)
%
% Coeffs:
%	A (k x k x q) array or a (k x k*q) matrix specifying the coefficient matrix
%	at each of the lags specified in the 'Lags' argument above. For a univariate
%	model, Coeffs can be a vector of length q.
%
% Lags:
%	Positive integer vector of the lags that correspond to each coefficient.
%	If this argument is not specified, it is taken to be the sequence (1:p),
%	where p is the number of slices in Coeffs.
%
% varargin:
% Name-Value pairs:
%	Name: 'J', Value: (Matrix) Coefficient of Y(t), [default: identity matrix]
%
% NOTES:
%	Eigensystem:
%	The Eigensystem of the Lag operator refers to the eigensystem of A, where
%	A is the matrix coefficient in the equivalent VAR(1) model for Y(t) = Phi(L)(Y(t)).
%
%	Transfer Function:
%	The Transfer Function, refers to the Transfer Function for the system:
%				w(t) = Y(t) - Phi(L)(Y(t))
%
% The following arithmetic operations can be performed on lag operators:-
%	1. (+, -) Unary and binary plus and minus.
%	2. (*) Multiplication.
%	3. (^) Power.
%	4. (==, ~= ), Elementwise comparison for equality (same as eq(), ne()).
%
% The other overloaded operations are:
%	5. (diag) Two lag operators, with matrix sizes k1 and k2 respectively,
%		can be combined along the diagonal to get an operator with size k1+k2.
%	6. plot
%	7. isequal
%	8. end
%	9. display
%
% Some of the other methods for LagOperator are:-
%	10. (Sub) Sub operators can be extracted.
%	11. (Apply, ()) The operator can be applied to time-series data.
%	12. (Evalues, Esystem, IsStable) The eigensystem can be computed and
%		stability can be determined.
%	13. (Permute) method permutes the operator to correspond to a permutation
%		in the order of the variables.
%	14. (Shift) The non-zero lags can be shifted by a scalar amount.
%	15. (TransferFunction) Calculates the associated transfer function.
%	16. (VAROne) Converts the operator to its equivalent AR1 form.
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
classdef (ConstructOnLoad) LagOperator
	
	properties (Access = 'private')
		Phi___ = [];
		Jmat__;
		lags = [];
	end
	
	properties (SetAccess = 'private')
		has_lag_zero = false;
	end
	
	properties (Access = 'private', Transient = true)
		% Per LagOperator object cache.
		hCache = [];
		b_store_as_compressed = false;
	end
	
	properties (Access = 'private', Dependent = true)
		Phi;
		Jmat;
	end
	

	methods (Static, Access = 'private')
		[lags, Coeffs, J, has_lag_zero] = ParseArgs(lags, Coeffs, J);
		[Coeffs, lags] = Trim(Coeffs, lags)
		M = BlockCompanion(C)
		M = BlockCompanionReduce(M, k)
		M = BlockCompanionSparse(C)
	end
	
	methods (Static)
		function obj = loadobj(obj)
			if isstruct(obj)
				if ~isfield(obj, 'Jmat__')
					obj.Jmat__ = [];
				end
				obj = LagOperator(obj.lags, obj.Phi___, 'J', obj.Jmat__);
			end
		end
		
		Obj = Concatenate(dim, Ops);
	end
	
	methods (Access = 'private')
		[left, right] = CheckBinaryOp__(left, right)
		Op = CheckUnaryOp__(Op)
		Op = Compose__(leftOp, rightOp)
	end
	
	methods
		% Usage:
		%	This = LagOperator(Lags, Coeffs, varargin)
		%
		% varargin:
		% Name-Value pairs:
		%	Name: 'J';			Value: (Matrix) Coefficient of Y(t)
		%
		function Op = LagOperator(Lags, Coeffs, varargin)
			if nargin == 0
				% Only when loaded from file.
				Op.hCache = Cache();
				return;
			end
			
			opts = parse_nvpairs(varargin, {'J', 'double', []});
			[Op.Phi, Op.lags, Op.Jmat, Op.has_lag_zero] = ParseArgs(Lags, Coeffs, opts.J);
			Op.hCache = Cache();
		end
		
		function J = get.Jmat(Op)
			J = Op.Jmat__;
			if issparse(J)
				J = full(J);
			end
		end
		
		function C = get.Phi(Op)
			C = Op.Phi___;
			if issparse(C)
				C = full(C);
				k = size(C, 1);
				q = size(C, 2);
				C = reshape(C, [k,k,q/k]);
			end
		end
		
		function Op = set.Jmat(Op, J)
			assert(ismatrix(J) && is_square(J));
			if Op.b_store_as_compressed
				Op.Jmat__ = sparse(J);
				
			elseif issparse(J)
				Op.Jmat__ = full(J);
				
			else
				Op.Jmat__ = J;
			end
		end
		
		function Op = set.Phi(Op, C)
			if Op.b_store_as_compressed
				Op.Phi___ = sparse(C(:,:));
				
			elseif issparse(C)
				Op.Phi___ = full(C);
				k = size(C, 1);
				q = size(C, 2);
				Op.Phi___ = reshape(Op.Phi__, [k,k,q/k]);
				
			else
				Op.Phi___ = C;
			end
		end
		
	end
		
	methods (Access = 'public')
		% Accessors
		Z = Apply(Op, Y, ivars_known_now, itime_range)
		a = CharacteristicPolynomial(This)
		[n, ijcounts] = CoeffCount(Op, eq_indices, varargin)
		[C, nzix, J] = Coeffs(Op, lags, varargin)
		[Q, L, Code] = Esystem(Op, Params)
		[l, ldistinct, lmulti] = Evalues(Op, Resolution)
		y = Generate(Op, N, y0, wn)
		n = InteractCount(Op, mask)
		bInteract = InteractMask(Op)
		b = IsEmpty(Op)
		[b, bInfinite, bHeuristic, loopexp, sr] = IsStable(Op, MyOne, bUseHeuristics, bSuppressDisplay)
		b = IsTrivial(Op)
		K = NumVars(This)
		[L, q] = Lags(Op)
		P = Order(Op)
		H = TransferFunction(Op, Freq)
		[sr, checkfailed] = SpectralRadius(Op)
		Op1 = VAROne(Op)
		
		% Overloaded Methods, all non-modifiers
		disp(Op)
		s = end(Op, k, n)
		iseq = eq(left, right)
		iseq = isequal(left, right)
		Opdiff = minus(left, right)
		Op = mpower(Op, n)
		OpProd = mtimes(left, right)
		isne = ne(left, right)
		g = plot(Oper, mLags, bInteracts, g, gParams, varargin)
		Opsum = plus(left, right)
		Opneg = uminus(Op)
		Op = uplus(Op)
		
		% Modifiers, must be called with Operator output argument.
		% Implementation note: Remember to clear or create new cache.
		Op = Denormalize(Op, lambda)
		[Op, lags] = Diag(Op)
		Op = Permute(Op, pseq)
		Op = Reduce(Op)
		opnew = Reflate(Op, sr_ratio)
		Op = Shift(Op, shft)
		[opnew, sr_old, sr_new] = Stabilize(Op, sr_new, sr_thresh, varargin)
		Op = StoreAsCompressed(Op, bstate);
		[Op, kseq] = Sub(Op, kseq, lagseq)
		
		% New Cache. Only 'This' object will get a new cache.
		function This = CacheNew(This)
			enforce_modifier(nargout);
			This.hCache = Cache();
		end
		
	end
	
end
%%
