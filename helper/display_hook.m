% display_hook
% Call this method to toggle display of any output of ObjType that occurs via DisplayFnString.
%
% Also see the Display_masker.m, disp_x.m
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
function display_hook(ObjType, DisplayFnString, X)
global truncation_length

% If overriding the builtins display and disp to call this function, then please
% note the following carefully:
%
% 1.
% Watch out when debugging! All debugging output to the command window appears
% to require calls to this method, so stepping through the code here will lead
% to a stack of function calls to this method, unless some means is found to
% turn off the output to the command window. Better to get @double/display.m right,
% copy here and replace 'double' by 'char'.
%
% 2.
% Found issue in openm.m
%	(\MATLAB\R2012a\toolbox\matlab\general\private\openm.m)
% There is no semicolon on line 12, after 'edit(filename)', causing 'display'
% to be called. This somehow leads to an error 'undefined command: display"
% whenever a file is opened for editing via the Current Folder window. Have
% no idea why. The function 'edit(filename)' just returns a boolean value
% so there is no reason to think that the missing semicolon is intentional.
% Hence edited the file openm.m (1/29/2013), in administrator mode, to add
% the absconding semicolon. This fixes the problem. 

	if Display_masker.GetState(ObjType) || Display_masker.GetState('all')
		
		% Suppress output.
		return;
		
	elseif defined('truncation_length')
		
		% Truncate output.
		[T, X] = evalc(X); % X is there to ensure that any expression gets evaluated.
		Ignore(X);
		n = min(truncation_length, length(T));
		builtin(DisplayFnString, T(1:n));
		
	else
		
		% Call builtin display.
		% Convert to temporary using ObjType to print result directly rather than doing X = ...
		
		if isa(X, 'function_handle') % function_handle(@xyz) doesn't work.
			disp(X);
			
		else
			htypefn = str2func(ObjType);
			builtin(DisplayFnString, htypefn(X));
		end
		
	end
	
end
%%
