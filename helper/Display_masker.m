% class Display_masker
% Suppresses output per object type using a persistent (hence shared) map associating object
% types with their suppression state (true or false). This should be used in conjuction with
% the disp_x() functions in .../Utility. These functions check the suppression
% map before calling the builtin functions display() and disp() respectively.
%
% Usage: First construct a Display_masker object:
%		dispmask = Display_masker(objtype, flg)
% If flg == true	: Upon construction, begins suppressing output of specified object type.
% If flg == false	: Upon construction, clears suppress flag for the specified object type.
% If flg is not specified, a Display_masker object is constructed without altering the
% the state of the suppression map.
%
% Methods:
% Flag(objtype, flg)
% Arguments:
%	objtype:
%		Builtin class name, builtin category or name of a user defined class that has a display
%		implementation. ex: 'char', 'double', 'MyClass', etc. To suppress output of all objects,
%		use 'all' as the object type.
%	flg:
%		logical true or false.
%
%	To start suppressing output for a particular object type:
%		dispmask.Flag(objtype, true)
%	This will globally suppress all output of the specified object type to the screen.
%
%	To stop suppressing output for a particular object type:
%		dispmask.Flag(objtype, false)
%	This will globally re-enable output of the specified object type to the screen.
%
%	To get the current supression state (true or false) of the specified object:
%		is_suppressed = dispmask.Flag(objtype);
%
% Returns:
%	The return value is the suppress state (logical true or false) of the object before any
%	change was made.
% __________________________________
%
%
% Restore()
%	Restores the previously set flag to its original state.
% __________________________________
%
%
% b = Display_masker.GetState(objtype)
%	Static method. Returns suppression state (true / false) of specified object type without
%	the need to construct a Display_masker object.
% __________________________________
%
% NOTES:
% 1. The Flag function can be called any number of times from within a function or script to
% suppress or re-enable output. For example, to suppress output from a particularly garrulous
% function while in a loop, call old = dispmask.Flag('all', true) before the chattering function
% and restore output after it using dispmask.Flags('all', old). A function that needs to
% ensure that its output is not suppressed can call Display_masker.Flags(objtype, false) as
% needed.
%
% 2. When a Display_masker object gets destroyed, it returns the suppression map to the state
% it was in prior to construction of that particular Display_masker object. So, assuming that
% the promises made by Matlab are true, the map will be returned to its original state at the
% end of its scope (which is function for matlab) regardless of the variations in the state
% within the scope. This cleanup will also occur in case of an error within the scope. This
% ensures that subsequent uses of Display_masker continue to work correctly.
%
% 3. Display_masker can be used for any object type provided the display() (or disp()) method
% for that type, or any other display wrapper or code, checks the flag using:
%		if Display_masker.GetState(objtype)
%			... no display ...
%		else
%			... display ...
%		end
%
% 4. Under certain circumstances (for example when you want to suppress output from functions
% not under your control) you may need to override the builtins, display() and disp() to include
% this check. However, this causes Matlab to trip up all over the place (regardless of how you
% do the override - short of cluttering all of your folders with private disp() functions). In
% view of this, for code under your control, the disp_x() option seems better and safer.
%
% 5. Using clear classes will reset the display mask to all false.
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
classdef Display_masker
	properties (Constant, Access = 'private')
		hmap = Structxh;
	end
	
	properties (Access = 'private')
		% See Matlab documentation on 'onCleanup'.
		cleanupObj;
		origdata = struct;
	end
	
	methods (Access = 'private')
		function Cleanup(This)
			This.Restore();
		end
	end
	
	methods (Static, Access = 'public')
		function currstate = GetState(objtype)
			[bval, bexists] = Display_masker.hmap.get_field(objtype);
			[bval_all, bexists_all] = Display_masker.hmap.get_field('all');
			currstate = (bexists && bval) || (bexists_all && bval_all);
		end
	end
	
	methods (Access = 'public')
		function This = Display_masker(objtypenm, flg)
			This.origdata = This.hmap.get_data();
			
			% See Matlab documentation on 'onCleanup'.
			This.cleanupObj = onCleanup(@() Cleanup(This)); % set _after_ setting origdata in 'This'.
			
			if (nargin == 2) && ~isempty(flg)
				This.Flag(objtypenm, flg);
			end
		end
		
		function Restore(This)
			This.hmap.set_data(This.origdata);
		end
		
		% Flag method should be non-static to ensure that a Display_masker object is constructed
		% before Flag() is called. The construction is required to ensure proper cleanup at end
		% of scope (scope is function for matlab) on error.
		function currstate = Flag(This, objtypenm, flg)
			assert( ischar(objtypenm) );
			currstate = Display_masker.GetState(objtypenm);
			if (nargin > 1)
				assert( islogical(flg) );
				if strcmpi(objtypenm, 'all') && (flg == false)
					This.hmap.clear_data();
				else
					This.hmap.set_field(objtypenm, flg);
				end
			end
		end

	end
	
end
%%
