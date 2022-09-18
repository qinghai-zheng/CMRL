function [sol,info] = prox_l12(x, gamma , param)
%PROX_L12 Proximal operator with L12 norm
%   Usage:  sol=prox_l12(x, gamma, param)
%           [sol,info] = prox_l12(x, gamma, param)
%
%   Input parameters:
%         x     : Input signal.
%         gamma : Regularization parameter.
%         param : Structure of parameters.
%   Output parameters:
%         sol   : Solution.
%         info  : Structure summarizing informations at convergence
%
%   prox_L12(x, gamma, param) solves:
%
%      sol = argmin_{z} 0.5*||x - z||_2^2 + gamma * ||x||_12^2
%
%   where 
%
%       ' || x ||_12 =  sum_j ( sum_i |x(i,j)|)^2  
%
%   The easiest way to use this proximla operator is to give a matrix x as
%   imput. In this case, the one norm will be computed over the lines (2nd
%   dimention) and the two norm will be computed over the rows (1st
%   dimention).
%
%   param is a Matlab structure containing the following fields: 
%
%    param.weights : weights for a weighted L12 norm (default = 1)
%
%    param.g_d, param.g_t are the group vectors. If you give a matrix,
%     do not set those parameters.
%
%     param.g_d contains the indices of the elements to be grouped and
%     param.g_t the size of the different groups. 
%
%     Warning: param.g_d and param.g_t have to be row vector!     
%     
%     Example: suppose x=[x1 x2 x3 x4 x5 x6] 
%                  and Group 1: [x1 x2 x4 x5] 
%                      group 2: [x3 x6]
%              
%     In matlab: 
%
%           param.g_d = [1 2 4 5 3 6]; param.g_t=[4 2];
%
%     Also this is also possible:
%
%           param.g_d = [4 5 3 6 1 2]; param.g_t=[2 4]; 
%
%    param.multi_group*: in order to group component in a not disjoint
%     manner, it is possible to use the multi_group option.
%     param.multi_group is now set automatically by the function. 
%
%     Overlaping group:
%     In order to make overlapping group just give a vector of g_d, g_b
%     and g_t. Example:
%       
%           param.g_d=[g_d1; g_d2; ...; g_dn];
%           param.g_t=[g_t1; g_t2; ...; g_tn];
%
%     Warning! There must be no overlap in g_d1, g_d2,... g_dn*
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
%   See also:  prox_l1 prox_linf1 prox_l21 prox_sumg
%
%   Demos: demo_compress_sensing3 demo_overlaping_groups_structured_sparsity
%
%   References:
%     F. Bach, R. Jenatton, J. Mairal, and G. Obozinski. Optimization with
%     sparsity-inducing penalties. arXiv preprint arXiv:1108.0775, 2011.
%     
%     M. Kowalski, K. Siedenburg, and M. Dorfler. Social sparsity!
%     neighborhood systems enrich structured shrinkage operators. Signal
%     Processing, IEEE Transactions on, 61(10):2498-2511, 2013.
%     
%     M. Kowalski. Sparse regression using mixed norms. Applied and
%     Computational Harmonic Analysis, 27(3):303-324, 2009.
%     
%     M. Kowalski and B. Torresani. Sparsity and persistence: mixed norms
%     provide simple signal models with dependent coefficients. Signal, image
%     and video processing, 3(3):251-264, 2009.
%     
%
%
%   Url: http://lts2research.epfl.ch/unlocbox/doc/prox/prox_l12.php

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
% Date: Mai 2013
% Testing: test_mixed_sparsity, test_tv

% Start the time counter
t1 = tic;

% Reshape x if not a row vector
t=size(x);

% Handle the parameter p
% argmin_x 1/2 ||x-y||_2^2 + p/gamma || x ||_{qp}^p
gamma = 2 * gamma;


% Optional input arguments
if nargin<3, param=struct; end

if ~isfield(param, 'g_d'),    param.g_d = 1:numel(x); end
if ~isfield(param, 'g_t')
    if numel(x) == size(x,1)*size(x,2); % matrix case
        param.g_t = size(x,2)*ones(1,size(x,1)); 
    else
        param.g_t = ones(1,numel(x)); 
    end
end

if ~isfield(param, 'weights'), param.weights=ones(numel(x),1) ; end
if ~isfield(param, 'verbose'), param.verbose=1 ; end
if ~isfield(param, 'multi_group'), param.multi_group=0 ; end



% Test of gamma
if test_gamma(gamma)
    sol = x;
    info.algo=mfilename;
    info.iter=0;
    info.final_eval=0;
    info.crit='--';
    info.time=toc(t1);
    return; 
end


% Test of the weights
param.weights=test_weights(param.weights);

test_multigroup(x,param.g_d,param.g_t);

% test if there are more than one group
if size(param.g_d,1)>1
    param.multi_group=1;
end



if param.multi_group==0
    
    % Number of group
    l=length(param.g_t);
 
   
    % Test if all the group have the same size
    if max(param.g_t)==min(param.g_t),
                
        Nelg = param.g_t(1);
        
        % reshape x in a useful manner
        X=transpose(x);
        X=X(param.g_d);
       
        X=reshape(X,numel(x)/l,l);
        % TODO: check transpose of the next line
        W=reshape(transpose(param.weights(param.g_d)),...
            numel(x)/l,l);
        
        [ Kw, ny] = find_mg( X,gamma,W );
        
        tau = gamma./( 1 + repmat(Kw',Nelg,1) * gamma) .* repmat(ny',Nelg,1);
        
        sol = soft_threshold( X, tau);
        
        
        % handle size for row vector
        if size(x,1)*size(x,2) == length(x)
            sol(param.g_d) = sol;
        else
            
            %reconstruct the solution
            sol(param.g_d)=transpose(sol);

        end

        
        
        
        
    else % group of different size
        x = x(:);
        sol=zeros(size(x));
        indice=0;
        
        xp=x;
        xp=xp(param.g_d);
        W=param.weights;
        W=W(param.g_d);
        for i=1:l
           temp=xp(indice+1:indice+param.g_t(i));
           w=W(indice+1:indice+param.g_t(i));
           [ Kw, ny] = find_mg( temp,gamma,w );
           tau = gamma./( 1 + gamma*Kw) * ny;
           s=soft_threshold( temp,tau);
           sol(indice+1:indice+param.g_t(i))=s;
           indice=indice+param.g_t(i);
           
        end
        %reconstruct the solution
        sol(param.g_d)=sol;
    end
    
    
    
    norm_L12=norm_l12(x,param.g_d,param.g_t,param.weights);
    % Log after the calculous of the prox
    if param.verbose >= 1
        fprintf('  prox_L12: ||x||_12^2 = %e\n', norm_L12^2);
    end

else % overlapping group
    r=size(param.g_t,1);
    
    % Parameter for the prox

    G=cell(r,1);
    
    for k=1:r
        param3.g_t=param.g_t(k,:);
        param3.g_d=param.g_d(k,:);
        param3.multi_group=0;
        param3.verbose=0;
        g.prox=@(x,T) prox_l12(x,T,param3);
        g.eval=@(x) norm_l12(x,param.g_d(k,:),param.g_t(k,:));
        G{k}=g;
    end
    
    param4.G=G;
    param4.verbose=param.verbose-1;
    sol=prox_sumg(x,gamma,param4);
    
end

%resahpe the solution
sol=reshape(sol,t);

iter=0;
crit='--';
info.algo=mfilename;
info.iter=iter;
info.final_eval=norm_l12(sol,param.g_d,param.g_t);
info.crit=crit;
info.time=toc(t1);

end

