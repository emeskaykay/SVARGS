%
%
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
function data = generate_data(mod, len)
	k = mod.opsparse.NumVars();
	p = mod.opsparse.Order();
	lead_count = max(25, 5*p);
	
	noise = mvnrnd(zeros(1, k), mod.sigma, (len + lead_count));
	noise = transpose(noise);
	noise = reshape(noise, [k, len + lead_count]);
	y_init = zeros(k, p);
	
	% Generate data values __________________
	data = mod.opsparse.Generate(lead_count + len, y_init, noise);
	data = data(:, lead_count+1:end);
	% _______________________________________
	
	% Adjust for any small mean remnant.
	data = bsxfun(@minus, data(:,:), mean(data(:,:), 2));
end

%%
