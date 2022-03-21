% Progress_display
% Display's progress with numeric value or string, without leaving a trail.
%
% CONSTRUCTOR ____________
% this = Progress_display(prefix, value_prefixes, varargin)
% prefix:
%	Unchanging string prefixes that identifies the task.
%
% value_prefixes:
%	String or cell-array of strings: {'<value name 1>', '<value name 2>', ... }
%	The task name will be printed before the updated numeric value or string.
%	The value names and values will be printed after the task name in the form:
%		<value name 1> = value_1, <value name 2> = value_2, ...
%	This argument should only be provided in the constructor. Subsequently
%	call the update() function with only the updated value or string.
%	Note that if no value names are provided, then calling update(vals) with
%	multiple values, will result in only the first value being printed.
%
% varargin:
% Optional name-value pairs
% name: 'value_formats'
% value(s): string or a cell array of strings that specifies the format for
%	each of the values in the prefix string. The format specification is the
%	same as in the function fprintf. This is only for formatting purposes.
%
% name: 'min_update_interval'
% value: Minimum update interval [default: 0.2 seconds]
%
% FUNCTIONS _____________
% function update_string(this, str)
%	Update the display string.
%
% update(this, val)
%	Update the numeric value displayed.
%
% increment(this, dn)
%	Increment the numeric value displayed.
%
% complete(this, final)
%	Call complete() after the concerned task is complete with the
%	final value to be displayed.
%
% NOTE 1:
%	If updates are made in very quick succession (say within 5 milliseconds)
%	the java text properties may not get updated in time. So for efficiency,
%	the update() method only updates the window if a minimum amount of time
%	has elapsed since the previous update. The minimum update interval is
%	a property (min_update_interval) that can be specified in the constructor.
%	Note that the last few updates may be skipped due to this behaviour.
%	Calling the complete() method at the end of the concerned task resolves
%	the skipped updates issue.
%
% NOTE 2:
%	To suppress display from Progress_display, the Display_masker('all') flag
%	must be set. See Display_masker.m.
%
% NOTE 3: Matlab handle class objects get _copied_ to and from the workers.
% So a Progress_display object in a parfor loop is deep copied to each of the
% workers.
%
% Backspace Idea: Thanks to Yair Altman's blog - UndocumentedMatlab.com
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
classdef Progress_display < handle

	properties (Constant)
		% Use handle object Strucxh to simulate a static non-constant variable.
		nesting = Structxh(struct('current_level', 0));
	end
	
	properties (Access = 'private')
		% CONSTANT
		min_update_interval = 0.2; % secs.
		last_line_length;
		task_name;
		value_prefixes;
		value_formats;
		commas;
		eqstrs;
		% _________
		
		% Handle to command window text area.
		jTextArea;
		
		curr_val;
		tlast;
		ambient_xdisplay_flag;
		orig_cmdedit_state;
		display_masker;
		inworker;
		in_use;
		is_queued;
		bforce;
		
		cleanupObj;
	end
	
	methods (Static, Access = 'private')
		% Uses Java methods.
		function jTextArea = GetCmdWndTextArea()
			jDesktop = com.mathworks.mde.desk.MLDesktop.getInstance;
			cmdWnd = jDesktop.getClient('Command Window');
			jTextArea = [];
			if ~isempty(cmdWnd)
				jTextArea = cmdWnd.getComponent(0).getViewport.getComponent(0);
			end
		end
		
		% Uses Java methods.
		function prev_edit_state = SetCmdWndKbEditable(jTextArea, beditable)
			if isempty(jTextArea)
				jTextArea = Progress_display.GetCmdWndTextArea();
			end
			
			prev_edit_state = [];
			if ~isempty(jTextArea) % Can happen if there is no window - for example in a parallel pool.
				prev_edit_state = jTextArea.isEditable();
				jTextArea.setEditable(beditable);
			end
		end
		
		% Uses Java methods.
		function line = GetCurrLine__(jTextArea)
			if isempty(jTextArea)
				jTextArea = Progress_display.GetCmdWndTextArea();
			end
			
			line = '';
			if ~isempty(jTextArea) % Can happen if there is no window - for example in a parallel pool.
				fulltxt = jTextArea.getText(); % Java String object.
				ixnewline = fulltxt.lastIndexOf(newline); % Returns -1 if not found.
				line = char(fulltxt.substring(ixnewline+1)); % No need to pass the end index.
			end
		end
		
		% Uses Java methods.
		function len = GetCurrLineLength__(jTextArea)
			line = Progress_display.GetCurrLine__(jTextArea);
			len = length(line);
		end
		
		function lev = increment_nesting_level__()
			lev = Progress_display.nesting.get_field('current_level');
			lev = lev + 1;
			Progress_display.nesting.set_field('current_level', lev);
		end
		
		function lev = decrement_nesting_level__()
			lev = Progress_display.nesting.get_field('current_level');
			lev = lev - 1;
			Progress_display.nesting.set_field('current_level', lev);
		end
		
		function k = ndigits__(n)
			k = 1 + ceil(log10(n));
		end
		
		function s = backspaces__(n)
			back = sprintf('\b');
			s = back(ones(1, n));
		end
		
		function cellstr = make_cell(strings)
			if ischar(strings)
				n = numel(strings);
				cellstr = cell(1, n);
				[cellstr{:}] = deal(strings);
				
			else
				assert(iscellstr(strings));
				cellstr = strings;
			end
			cellstr = vecr(cellstr);
		end
		
	end
	
	methods (Access = 'private')
		function Cleanup(this)
			Progress_display.decrement_nesting_level__();
			if ~ishandle(this)
				Progress_display.SetCmdWndKbEditable([], true);
			else
				Progress_display.SetCmdWndKbEditable(this.jTextArea, this.orig_cmdedit_state);
			end
		end
		
		function reset__(this)
			this.curr_val = 0;
			if ~isempty(this.orig_cmdedit_state)
				Progress_display.SetCmdWndKbEditable(this.jTextArea, this.orig_cmdedit_state);
			end
			this.jTextArea = Progress_display.GetCmdWndTextArea();
			this.reset_tlast__();
			this.inworker = [];
			this.in_use = 0;
		end
		
		function b = is_inside_worker__(this)
			if isempty(this.inworker)
				this.inworker = is_worker();
			end
			b = this.inworker;
		end
		
		function b = suppress_display__(this)
			if this.is_inside_worker__() || ((this.in_use > 1) && ~this.is_queued)
				b = true;
				
			elseif this.bforce
				b = this.ambient_xdisplay_flag;
				
			else
				b = (this.ambient_xdisplay_flag || (etime(clock(), this.tlast) < this.min_update_interval));
			end
		end
		
		function valstr = make_value_string__(this, vals)
			n = min(numel(vals), numel(this.value_prefixes));
			
			if isnumeric(vals)
				str = [this.value_prefixes(1:n); this.eqstrs(1:n); this.value_formats(1:n); this.commas(1:n)];
				fmtstr = [str{:}];
				fmtstr = fmtstr(1:end-2);
				valstr = num2str(vals, fmtstr);
			else
				assert(iscell(vals));
				valstr = [this.value_prefixes(1:n); this.eqstrs(1:n); vals(1:n); this.commas(1:n)];
				valstr = valstr(1:end-2);
			end
		end
		
		function do_update__(this, vals)
			valstr = make_value_string__(this, vals);
			this.do_update_string__(valstr);
		end
		
		function line_length = do_update_string__(this, str)
			assert(ischar(str));
			
			% Temporarily set mask to ambient display state.
			this.display_masker.Restore();
			
			dispstr = [this.task_name, ': ', str];
			line_length = Progress_display.GetCurrLineLength__(this.jTextArea);
			
			if (line_length == 0)
				print_x('%s', dispstr);
				this.last_line_length = length(dispstr);
				
			elseif isempty(this.last_line_length) || (line_length == this.last_line_length)
				print_x('%s', [this.backspaces__(this.last_line_length), dispstr]);
				this.last_line_length = length(dispstr);
				
			else
 				this.nextline(' Progress display interrupted ... ');
 				print_x('%s', dispstr);
				this.last_line_length = length(dispstr);
			end
			
			% Re-suppress all calls to display.
			this.display_masker.Flag('all', true);
			
			this.tlast = clock();
		end
		
		function reset_tlast__(this)
			this.tlast = clock();
			this.tlast(1) = 1;
		end
	end
	
	methods
		function this = Progress_display(prefix, value_prefixes, varargin)
			this.cleanupObj = onCleanup(@() Cleanup(this));
			
			assert(ischar(prefix));
			if (nargin < 2) || isempty(value_prefixes), value_prefixes = ''; end
			if ischar(value_prefixes), value_prefixes = {value_prefixes}; end
			assert(iscellstr(value_prefixes));
			
			reset__(this);
			
			this.task_name = prefix;
			this.value_prefixes = value_prefixes;
			opts = parse_nvpairs(varargin, ...
				{'value_formats', '', '%g'}, ...
				{'min_update_interval', 'double', 0.2}, ...
				{'is_queued', 'logical', false});
			
			assert(isscalar(opts.min_update_interval) && (opts.min_update_interval > 0));
			this.min_update_interval = opts.min_update_interval;
			this.is_queued = opts.is_queued;
			
			n = numel(this.value_prefixes);
			if ischar(opts.value_formats)
				valfmt = opts.value_formats;
				opts.value_formats = cell(1, n);
				[opts.value_formats{:}] = deal(valfmt);
			end
			this.value_formats = opts.value_formats;
			assert(iscell(this.value_formats) && (numel(this.value_formats) == n));
			
			this.commas = cell(1, n);
			[this.commas{:}] = deal(', ');
			
			this.eqstrs = cell(1, n);
			[this.eqstrs{:}] = deal(' = ');
			
			% Save ambient display state and Suppress all calls to xdisplay.
			this.display_masker = Display_masker();
			this.ambient_xdisplay_flag = this.display_masker.Flag('all', true);
			
			this.bforce = false;
			this.in_use = Progress_display.increment_nesting_level__();
			
			this.orig_cmdedit_state = Progress_display.SetCmdWndKbEditable(this.jTextArea, false);
		end
		
		function update_string(this, str, bforce_once)
			if (nargin < 3) || isempty(bforce_once), bforce_once = false; end
			
			if ischar(str), str = {str};
			else
				assert(iscell(str));
			end
			
			this.bforce = bforce_once;
			if ~this.suppress_display__()
				str = make_value_string__(this, str);
				this.do_update_string__(str);
			end
			this.bforce = false;
		end
		
		function update(this, val, bforce_once)
			if (nargin < 3) || isempty(bforce_once), bforce_once = false; end
			assert(isnumeric(val));
			
			this.bforce = bforce_once;
			this.curr_val = val;
			if ~this.suppress_display__()
				this.do_update__(val);
			end
			this.bforce = false;
		end
		
		function increment(this, dn)
			assert(isnumeric(dn));
			
			this.curr_val = this.curr_val + dn;
			if ~this.suppress_display__()
				this.do_update__(this.curr_val);
			end
		end
		
		function force_next(this)
			this.bforce = true;
		end
		
		function nextline(this, msgstr)
			if (nargin < 2), msgstr = ''; end
			assert(ischar(msgstr));
			
			% Temporarily set mask to ambient display state.
			this.display_masker.Restore();
			
			print_x('%s', [msgstr, newline]);
			
			% Re-suppress all calls to xdisplay.
			this.display_masker.Flag('all', true);
			
			this.reset_tlast__();
		end
		
		function complete(this, final_val)
			if (nargin < 2) || isempty(final_val)
				% do nothing
			elseif isnumeric(final_val)
				this.curr_val = final_val;
				this.update(final_val, true);
			else
				this.update_string(final_val, true);
			end
			
			this.display_masker.Restore();
			muisaved = this.min_update_interval;
			this.min_update_interval = 0;
			if ~this.suppress_display__()
				print_x(' %s\n', [this.task_name, ' Completed. ']);
			end
			this.min_update_interval = muisaved;
			
			this.reset__();
			this.cleanupObj = [];
		end
	end
	
end
%%

