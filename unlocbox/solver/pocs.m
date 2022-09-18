function [sol, info,objective] = pocs(x_0,F, param)
%POCS Projection onto convex sets
%   Usage: sol = pocs(x_0,F, param);
%          sol = pocs(x_0,F);
%          [sol,info,objective] = pocs(...);
%
%   Input parameters:
%         x_0   : Starting point of the algorithm
%         F     : Array of function to minimize
%         param : Optional parameter
%   Output parameters:
%         sol   : Solution
%         info  : Structure summarizing informations at convergence
%         objective: vector (evaluation of the objective function each iteration)
%
%   backard_backward solves:
%
%      sol = argmin || x - x_0 ||_2     for x belong to R^N and belonging to all sets
%
%   where x is the variable.
%
%    x_0 is the starting point.
%
%    F is an array of structures representing functions containing 
%     operators inside and eventually the norm. The prox: F(i).prox and 
%     the function: F(i).eval are defined in the same way as in the 
%     Forward-backward and Douglas-Rachford algorithms F(i).prox should 
%     perform the projection. This .prox notation is kept for compatibility
%     reason.
%
%    param a Matlab structure containing the following fields:
%
%     General parameters:
%
%      param.tol : is stop criterion for the loop. The algorithm stops if
%
%           (  n(t) - n(t-1) )  / n(t) < tol,
%      
%       where  n(t) = f_1(x)+f_2(x) is the objective function at iteration t*
%       by default, tol=10e-4.
%
%      param.maxit : is the maximum number of iteration. By default, it is 200.
% 
%      param.verbose : 0 no log, 1 print main steps, 2 print all steps.
%
%      param.abs_tol : If activated, this stopping critterion is the
%       objectiv function smaller than param.tol. By default 1.
%
%
%   info is a Matlab structure containing the following fields:
%
%    info.algo : Algorithm used
%
%    info.iter : Number of iteration
%
%    info.time : Time of exectution of the function in sec.
%
%    info.final_eval : Final evaluation of the objectivs functions
%
%    info.crit : Stopping critterion used 
%
%    info.rel_norm : Relative norm at convergence 
%
%
%   See also:  backward_backward generalized_forward_backward
%
%   Url: http://lts2research.epfl.ch/unlocbox/doc/solver/pocs.php

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



% Author: Nathanael Perraudin
% Date: 14 dec 2012

% Start the time counter
t1 = tic;

% test the evaluate function
[F] = test_eval(F);

% number of function
m = size(F,2);

% Optional input arguments
if nargin<3, param=struct; end

if ~isfield(param, 'tol'), param.tol=10e-4 ; end
if ~isfield(param, 'maxit'), param.maxit=200; end
if ~isfield(param, 'verbose'), param.verbose=1 ; end
if ~isfield(param, 'abs_tol'), param.abs_tol=1 ; end





% Initialization    
curr_norm = 0;
for ii = 1:m
    curr_norm = F{ii}.eval(x_0)+curr_norm;
end
[~,~,prev_norm,iter,objective,~] = convergence_test(curr_norm);

sol=x_0;

% Main loop
while 1
    
    %
    if param.verbose >= 2
        fprintf('Iteration %i:\n', iter);
    end
    
    

    for ii = 1:m
       sol = F{ii}.prox(sol,0);
    end

    
    % Global stopping criterion
    curr_norm=0;
    for ii = 1:m
       curr_norm = F{ii}.eval(sol)+curr_norm;
    end
    [stop,rel_norm,prev_norm,iter,objective,crit] = convergence_test(curr_norm,prev_norm,iter,objective,param);
    [sol,param]=post_process(sol,iter,curr_norm,prev_norm,param);
    if stop
        break;
    end
    if param.verbose >= 2
        fprintf('  ||f|| = %e, rel_norm = %e\n', ...
            curr_norm, rel_norm);
    end

    
end

% Log
if param.verbose>1
    % Print norm
    fprintf('\n Solution found:\n');
    fprintf(' Final residue:%e     Final relative norm: %e\n',curr_norm, rel_norm );
    
    
    % Stopping criterion
    fprintf(' %i iterations\n', iter);
    fprintf(' Stopping criterion: %s \n\n', crit);
    
end

if param.verbose==1
    % single line print
    fprintf('  POCS: ||f||=%e, REL_OBJ=%e, iter=%i, crit: %s \n',curr_norm, rel_norm,iter,crit );
end


info.algo=mfilename;
info.iter=iter;
info.final_eval=curr_norm;
info.crit=crit;
info.time=toc(t1);
info.rel_norm=rel_norm;

end
