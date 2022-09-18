function [sol,info] = prox_l2gradfourier(x, gamma, param)
%PROX_L2gradfourier Proximal operator of the 2 norm of the gradient in the Fourier domain
%   Usage:  sol=prox_l2gradfourier(x, gamma)
%           sol=prox_l2gradfourier(x, gamma, param)
%           [sol, info]=prox_l2gradfourier(x, gamma, param)
%
%   Input parameters:
%         x     : Input signal.
%         gamma : Regularization parameter.
%         param : Structure of optional parameters.
%   Output parameters:
%         sol   : Solution.
%         info  : Structure summarizing informations at convergence
%
%   This function compute the 1 dimensional proximal operator of x. For
%   matrices, the function is applied to each column. The parameter
%   param.d2 allows the user to use the 2 dimentional gradient.
%
%   Warning: the signal should not be centered. Indice 1 for abscissa 0.
%
%   PROX_L2GRADFOURIER(x, gamma, param) solves:
%
%      sol = argmin_{z} 0.5*||x - z||_2^2 + gamma * ||grad(Fz)||_2^2
%
%   param is a Matlab structure containing the following fields:
%
%    param.weights : weights if you use a an array.
%
%    param.verbose : 0 no log, 1 a summary at convergence, 2 print main
%     steps (default: 1)
%
%    param.deriveorder : Order ot the derivative default 1
%
%    param.d2 : 2 dimentional gradient (default 0)
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
%   See also:  prox_l2 prox_l2grad prox_tv
%
%   Url: http://lts2research.epfl.ch/unlocbox/doc/prox/prox_l2gradfourier.php

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
% Date: Jan 2013
%

% Start the time counter
t1 = tic;

% Optional input arguments
if nargin<3, param=struct; end

% Optional input arguments
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'weights'), param.weights = 1; end
if ~isfield(param, 'deriveorder'), param.deriveorder = 1; end
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


p=param.deriveorder;

% useful function
h=@(t) 1./(1+param.weights(:)*gamma*t).^p;

if param.d2
    L=size(x,1);
    Q=size(x,2);
    l=(0:L-1)';
    q=(0:Q-1);
    lambda=(2-2*cos(2*pi*l/L))*(2-2*cos(2*pi*q/Q));
    sol=x.*repmat(h(lambda),[1,1,size(x,3)]);
    
    [dx, dy] =  gradient_op(1/sqrt(L)*1/sqrt(Q)*fft2(sol));
    curr_norm = norm(dx.^2+dy.^2,'fro')^2;
elseif (size(x,1)==1) || (size(x,2)==1)
    L=size(x,1)*size(x,2);
    if size(x,1)>size(x,2)
        l=(0:L-1)';
    else
        l=(0:L-1);
    end
    lambda=2-2*cos(2*pi*l/L);        

    %filtering
    sol=x.*h(lambda);
    curr_norm = norm(gradient_op1d(1/sqrt(L)*fft(sol(:))),'fro')^2;
else
    L=size(x,1);
    l=(0:L-1)';
    lambda=2-2*cos(2*pi*l/L);
    sol = x.*repmat(h(lambda),1,size(x,2));
    curr_norm = norm(gradient_op1d(1/sqrt(L)*fft(sol)),'fro')^2;
end

% Summary
if param.verbose>=1
   fprintf('  Prox_l2gradfourier: 1 iteration, ||grad(Fx)||^2=%g\n',curr_norm);
end

% zero iteration
iter=1;
crit='--';
info.algo=mfilename;
info.iter=iter;
info.final_eval=curr_norm;
info.crit=crit;
info.time=toc(t1);

end




