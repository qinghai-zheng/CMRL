%DEMO_ADMM  Example of use of the admm solver
%
%   We present an example of the admm solver through an image
%   reconstruction problem.
%   The problem can be expressed as this
%
%       argmin ||x||_TV s.t ||b-Ax||_2 < epsilon
%
%   Note that the constraint can be inserted in the objective function
%   thanks to the help of the indicative function. Then we recover the
%   general formulation used for the solver of this toolbox.
%
%   In order to solve this problem with the admm, we need to transform 
%   the problem a little. So, we formulate it in a different way.
%
%       argmin ||x||_TV s.t ||b-Ax||_2 < epsilon
%
%   Where b is the degraded image, I the identity and A an operator representing the mask. We set 
%
%    f_1(x)=||x||_{TV}
%     We define the prox of f_1 as: 
%
%        prox_{f1,gamma} (z) = argmin_{x} 1/2 ||x-z||_2^2  +  gamma ||z||_TV
%
%    f_2 is the indicator function of the set S define by Ax-b||_2 < epsilon
%     We define the prox of f_2 as 
%
%        prox_{f2,gamma} (z) = argmin_{x} 1/2 ||x-z||_2^2  +  gamma i_S( x ),
%
%     with i_S(x) is zero if x is in the set S and infinity otherwise.
%     This previous problem has an identical solution as:
%
%        argmin_{z} ||x - z||_2^2   s.t.  ||b - A z||_2 < epsilon
%
%     It is simply a projection on the B2-ball.
%
%   Results
%   -------
%
%   Figure 1: Original image
%
%      This figure shows the original Lena image. 
%
%   Figure 2: Depleted image
%
%      This figure shows the image after the application of the mask. Note
%      that 70% of the pixels have been removed.
%
%   Figure 3: Reconstructed image
%
%      This figure shows the reconstructed  image thanks to the algorithm.
%   
%
%   References:
%     P. Combettes and J. Pesquet. Proximal splitting methods in signal
%     processing. Fixed-Point Algorithms for Inverse Problems in Science and
%     Engineering, pages 185-212, 2011.
%     
%
%   Url: http://lts2research.epfl.ch/unlocbox/doc/demos/demo_admm.php

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
%


%% Initialisation

clear all;
close all;

% Loading toolbox
init_unlocbox();

verbose = 2;    % verbosity level

%% Defining the problem

tau = 1;        % regularization parameter

% Original image
im_original = lena();
   
% Creating the problem
A = rand(size(im_original));
A = (A > 0.7);

% Depleted image
b = A .* im_original;

%% Defining proximal operators
% Define the prox of f2 see the function proj_B2 for more help
operatorA = @(x) A .* x;
param_proj.epsilon = 0;
param_proj.A = operatorA;
param_proj.At = operatorA;
param_proj.y = b;
param_proj.verbose = verbose -1;

% setting the function f2 
f2.prox = @(x, T) proj_b2(x, T, param_proj);
f2.eval = @(x) norm(A(:) .* x(:) - b(:))^2;

% setting the function f1 (norm TV)
param_tv.verbose = verbose - 1;
param_tv.maxit = 50;

f1.prox = @(x, T) prox_tv(x, tau*T, param_tv);
f1.eval = @(x) tau * norm_tv(x);   

% setting the frame L
L = @(x) x;
% L= eye(size(b)); % Eventually, this choice is possible (Handle carefully)


%% solving the problem

% setting different parameter for the simulation
param_solver.verbose = verbose;     % display parameter
param_solver.maxit = 100;           % maximum iteration
param_solver.tol = 1e-3;            % tolerance to stop iterating
param_solver.gamma = 0.1;           % stepsize
                                    % Be careful with this parameter

param_solver.L = L;                 % Special paramter for ADMM

sol = admm(b, f1, f2, param_solver);

%% displaying the result
imagesc_gray(im_original, 1, 'Original image');
imagesc_gray(b, 2, 'Depleted image');
imagesc_gray(sol, 3, 'Reconstructed image');
    

%% Closing the toolbox
close_unlocbox();

