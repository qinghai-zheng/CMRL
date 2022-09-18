function [sol,info] = prox_l2grad(x, gamma, param)
%PROX_L2grad Proximal operator of the 2 norm of the gradient in 1 dimension
%   Usage:  sol=prox_l2grad(x, gamma)
%           sol=prox_l2grad(x, gamma, param)
%           [sol, info]=prox_l2grad(x, gamma, param)
%
%   Input parameters:
%         x     : Input signal.
%         gamma : Regularization parameter.
%         param : Structure of optional parameters.
%   Output parameters:
%         sol   : Solution.
%         infos : Structure summarizing informations at convergence
%
%   This function compute the 1 dimensional proximal operator of x. For
%   matrices, the function is applied to each column. To use the 2D
%   proximal operator just set up the parameter param.2d to 1.
%
%   PROX_L2GRAD(x, gamma, param) solves:
%
%      sol = argmin_{z} 0.5*||x - z||_2^2 + gamma * ||grad(Az)||_2^2
%
%   param is a Matlab structure containing the following fields:
%
%    param.abasis : to use another basis than the DFT (default: 0). To be
%                     done -- Not working yet
%
%    param.weights : weights if you use a an array.
%
%    param.verbose : 0 no log, 1 a summary at convergence, 2 print main
%     steps (default: 1)
%
%    param.d2 : 2 dimentional gradient (default 0)
%
%    param.A : Forward operator (default: Id).
%
%    param.At : Adjoint operator (default: Id).
%
%    param.tight : 1 if A is a tight frame or 0 if not (default = 1)
%
%    param.nu : bound on the norm of the operator A (default: 1), i.e.
%
%        ` ||A x||^2 <= nu * ||x||^2 
%
%    param.tol : is stop criterion for the loop. The algorithm stops if
%
%         (  n(t) - n(t-1) )  / n(t) < tol,
%      
%     where  n(t) = f(x)+ 0.5 X-Z_2^2 is the objective function at iteration t*
%     by default, tol=10e-4.
%
%    param.maxit : max. nb. of iterations (default: 200).
%
%    param.deriveorder : Order ot the derivative default 1
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
%    info.final_eval : Final evaluation of the function
%
%    info.crit : Stopping critterion used 
%
%
%   See also:  prox_l2 prox_tv prox_l2gradfourier
%
%   Url: http://lts2research.epfl.ch/unlocbox/doc/prox/prox_l2grad.php

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
% Date: Nov 2012
%

% Start the time counter
t1 = tic;

% Optional input arguments
if nargin<3, param=struct; end

% Optional input arguments
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'weights'), param.weights = 1; end
if ~isfield(param, 'deriveorder'), param.deriveorder = 1; end
if ~isfield(param, 'abasis'), param.abasis = 0; end
if ~isfield(param, 'tight'), param.tight = 1; end
if ~isfield(param, 'nu'), param.nu = 1; end
if ~isfield(param, 'tol'), param.tol = 1e-3; end
if ~isfield(param, 'maxit'), param.maxit = 200; end
if ~isfield(param, 'At'), param.At = @(x) x; end
if ~isfield(param, 'A'), param.A = @(x) x; end
if ~isfield(param, 'd2'), param.d2 = 0; end

warning=0;
% test the parameters
if test_gamma(gamma)
    sol = x;
    info.algo=mfilename;
    info.iter=0;
    info.final_eval=0;
    info.crit='--';
    info.time=toc(t1);
    return; 
end

test_weights(param.weights,warning);

if param.abasis
    error('This option is not currently supported, please contact the developer')
end

p=param.deriveorder;

if param.tight
    % useful function
    h=@(t) 1./(1+param.weights*param.nu*gamma*t).^p;

    % size of the signal

    temp=param.A(x);

    

    if param.d2
        
        L=size(temp,1);
        Q=size(temp,2);
        l=(0:L-1)';
        q=(0:Q-1);
        lambda=(2-2*cos(2*pi*l/L))*(2-2*cos(2*pi*q/Q));

        
        sol=x+1/param.nu*param.At(ifft2(fft2(temp).*...
            repmat(h(lambda),[1,1,size(temp,3)]))-temp);
                
        [dx, dy] =  gradient_op(param.A(sol)  );
        curr_norm = norm(dx.^2+dy.^2,'fro')^2;
        
    elseif (size(temp,1)==1) || (size(temp,2)==1)
        L=size(temp,1)*size(temp,2);
        if size(x,1)>size(x,2)
            l=(0:L-1)';
        else
            l=(0:L-1);
        end
        lambda=2-2*cos(2*pi*l/L);        

        %filtering
        sol=x+1/param.nu*param.At(ifft(fft(temp).*h(lambda))-temp);
        curr_norm = norm(gradient_op1d(reshape(param.A(sol),[],1)),'fro')^2;
    else
        L=size(x,1);
        dim = size(x,2);
        l=(0:L-1)';
        lambda=2-2*cos(2*pi*l/L);        

        %filtering
        sol=x+1/param.nu*param.At(ifft(fft(temp).* ...
            repmat(h(lambda),1,dim))-temp); 
        curr_norm = norm(gradient_op1d(param.A(sol)),'fro')^2;

     end



    % one iteration
    iter=1;
    crit='--';
else % non tight frame case (gradient descent)
     error('Not done yet!!! Sorry')
  
end

% Check if x is real
if isreal(x)
   sol = real(sol); 
end

% Summary
if param.verbose>=1
   fprintf('  Prox_l2grad: %i iteration(s), ||grad(x)||^2=%g\n',iter,curr_norm);
end



info.algo=mfilename;
info.iter=iter;
info.final_eval=curr_norm;
info.crit=crit;
info.time=toc(t1);

end




