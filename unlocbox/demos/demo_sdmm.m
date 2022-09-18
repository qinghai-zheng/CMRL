%DEMO_SDMM  Example of use of the sdmm solver 
%
%   We present an example of the solver through an image
%   denoising problem.
%   The problem can be expressed as this
%
%       argmin_x,y,z ||x-b||_2^2 + tau1*||y||_TV + tau2 * ||H(z)||_1 such that  x = y = Hz
%
%   Where b is the degraded image, tau_1 and tau_2 two real positive constant and H a linear operator on x.
%   H is a wavelet operator. We set:
%
%    g_1(x)=||x||_{TV}
%     We define the prox of g_1 as: 
%
%        prox_{f1,gamma} (z) = argmin_{x} 1/2 ||x-z||_2^2  +  gamma ||z||_TV
%
%    g_2(x)=||H(x)||_1
%     We define the prox of g_2 as:
%
%        prox_{f1,gamma} (z) = argmin_{x} 1/2 ||x-z||_2^2  +  gamma ||H(z)||_1
%
%    f(x)=||x-b||_2^2
%     We define the gradient as: 
%
%        grad_f(x) = 2 * (x-b)
%
%   Results
%   -------
%
%   Figure 1: Original image
%
%      This figure shows the original image (The cameraman). 
%
%   Figure 2: Depleted image
%
%      This figure shows the image after addition of the noise
%
%   Figure 3: Reconstruted image
%
%      This figure shows the reconstructed image thanks to the algorithm.
%   
%   The rwt toolbox is needed to run this demo.
%
%   References:
%     P. Combettes and J. Pesquet. Proximal splitting methods in signal
%     processing. Fixed-Point Algorithms for Inverse Problems in Science and
%     Engineering, pages 185-212, 2011.
%     
%
%   Url: http://lts2research.epfl.ch/unlocbox/doc/demos/demo_sdmm.php

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
% Date: November 2012

%% Initialisation

clear all;
close all;

% Loading toolbox
init_unlocbox();

verbose = 2;    % verbosity level

%% Defining the problem

% Original image
im_original = lena();

% Depleted image
sigma=0.9;
b=im_original+sigma^2*rand(size(im_original));

%% Defining proximal operators

% setting different parameter for the simulation
tau1 = 0.2; % regularization parameter for the TV norm
tau2 = 0.2; % regularization parameter for the wavelet

% setting the function f1 

% -----  TV norm ------
param_tv.verbose = verbose - 1;
param_tv.maxit = 100;

g_tv.prox = @(x, T) prox_tv(x, T*tau1, param_tv);
g_tv.eval = @(x) tau1 * norm_tv(x);   

g_tv.x0 = b;
g_tv.L = @(x) x;
g_tv.Lt = @(x) x;

% -----  Wavelet ------
L=8;
% h = daubcqf(2);
% A2 = @(x) mdwt(x,h,L);
% A2t = @(x) midwt(x,h,L);
A2 = @(x) fwt2(x,'db1',L);
A2t = @(x) ifwt2(x,'db1',L);

param_l1.verbose = verbose - 1;

g_l1.prox = @(x, T) prox_l1(x, T*tau2, param_l1);
g_l1.eval = @(x) tau2*norm(reshape(A2(x),[],1),1);  

g_l1.L = A2;
g_l1.Lt = A2t;
g_l1.x0 = b;


% -----  L2 norm ------
g_l2.prox=@(x, T) (2*T*b+x)/(2*T+1);
g_l2.eval=@(x) eps;  
g_l2.x0=b;
g_l2.L=@(x) x;
g_l2.Lt=@(x) x;

%% Solving the problem

% Parameter for the sum of function: F
F={g_l1,g_tv, g_l2};
param4.maxit=30;

% solving the problem
sol=sdmm(F,param4);


%% displaying the result
imagesc_gray(im_original, 1, 'Original image');
imagesc_gray(b, 2, 'Depleted image');
imagesc_gray(sol, 3, 'Reconstructed image');
    

%% Closing the toolbox
close_unlocbox();


