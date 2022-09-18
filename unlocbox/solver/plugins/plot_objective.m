function [ sol ] = plot_objective(x, fig,setlim)
%PLOT_OBJECTIVE Plot objective function over iters(plugin for the UNLocBoX)
%   Usage [ im ] = plot_image( im, fig );
%         [ im ] = plot_image( im, fig , setlim);
%
%   Input parameters:
%         info_iter   : Structure of info of current iter of algorithm
%         fig   : Figure
%         setlim : Set limite for the image (default 0)
%
%   Output parameters:
%         im    : Input image
%
%   This plugin display the image every iterations of an algorithm. To use
%   the plugin juste define:
%       
%       fig=figure(100);
%       param.do_sol=@(x) plot_image(x,fig);
%
%   In the structure of optional argument of the solver.
%
%
%   Url: http://lts2research.epfl.ch/unlocbox/doc/solver/plugins/plot_objective.php

% Copyright (C) 2012-2013 Nathanael Perraudin.
% This file is part of UNLOCBOX version 1.5.0
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% Author: Vassilis Kalofolias
% Date  : November 2013


if nargin<3
    setlim = 0; 
end

% select the figure
if x.iter<2
    figure(fig);
end

sol = x.sol;
sol(sol<0) = 0;
sol(sol>1) = 1;

% display the image
subplot 211
imshow(sol, []);
title(['Current it: ', num2str(x.iter),'   Curr obj: ', ...
    num2str(x.curr_norm)]);
subplot 212
semilogy(x.objective); title('Objective function')
drawnow;

% return the image
if setlim
   x.sol = sol; 
end
sol=x.sol;

end


