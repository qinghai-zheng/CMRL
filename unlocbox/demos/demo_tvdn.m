%DEMO_TVDN Demonstration of the use of the tvdn solver
%
%   In this demo we solve two different problems. Both can be written on this form:
%
%      sol arg min ||x||_TV  s.t.  ||y-A x||_2 < epsilon
%
%   The first problem is an inpainting problem with 33% of the pixel. In
%   that case A is simply a mask and y the know pixels.
%
%   The second problem consists of reconstructing the image with only 33%
%   of the Fourier coefficients. In that case A is a truncated Fourier
%   operator.
% 
%   Figure 1: Original image
%
%      The cameraman
%
%   Figure 2: First problem
%
%      In painting with 33% of known pixel and a SNR of 30dB
%
%   Figure 3: Second problem
%
%      Image reconstruction with 33% of the Fourier measurements and a SNR
%      of 30dB.
%
%   Url: http://lts2research.epfl.ch/unlocbox/doc/demos/demo_tvdn.php

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

% Author: Gilles Puy, Nathanael Perraudin
% Date: Nov. 1, 2012


%% Initialisation

clear all;
close all;

% Loading toolbox
init_unlocbox();

verbose = 2;    % verbosity level

%% Parameters
input_snr = 30; % Noise level (on the measurements)

%% Load an image
im = cameraman();
imagesc_gray(im, 1, 'Original image');


%% Create a mask with 33% of 1 (the rest is set to 0)

% Mask
mask = rand(size(im)) < 0.33; ind = find(mask==1);
% Masking matrix (sparse matrix in matlab)
Ma = sparse(1:numel(ind), ind, ones(numel(ind), 1), numel(ind), numel(im));
% Masking operator
A = @(x) Ma*x(:); % Select 33% of the values in x;
At = @(x) reshape(Ma'*x(:), size(im)); % Adjoint operator + reshape image

%% First problem: Inpainting problem
% Select 33% of pixels
y = A(im);
% Add Gaussian i.i.d. noise
sigma_noise = 10^(-input_snr/20)*std(im(:));
y = y + randn(size(y))*sigma_noise;
% Display the downsampled image
imagesc_gray(At(y),2,'Measured image',121);
% Parameters for TVDN
param.verbose = 1; % Print log or not
param.gamma = 1e-1; % Converge parameter
param.tol = 1e-4; % Stopping criterion for the TVDN problem
param.maxit = 200; % Max. number of iterations for the TVDN problem
param.maxit_tv = 100; % Max. nb. of iter. for the sub-problem (proximal TV operator)
param.nu_b2 = 1; % Bound on the norm of the operator A
param.tol_b2 = 1e-4; % Tolerance for the projection onto the L2-ball
param.tight_b2 = 0; % Indicate if A is a tight frame (1) or not (0)
param.maxit_b2 = 500;
% Tolerance on noise
epsilon = sqrt(chi2inv(0.99, numel(ind)))*sigma_noise;
% Solve TVDN problem
sol = solve_tvdn(y, epsilon, A, At, param);
% Show reconstructed image
imagesc_gray(sol, 2, 'Reconstructed image',122);


%% Second problem: Reconstruct from 33% of Fourier measurements
% Composition (Masking o Fourier)
A = @(x) Ma*reshape(fft2(x)/sqrt(numel(im)), numel(x), 1);
At = @(x) ifft2(reshape(Ma'*x(:), size(im))*sqrt(numel(im)));
% Select 33% of Fourier coefficients
y = A(im);
% Add Gaussian i.i.d. noise
sigma_noise = 10^(-input_snr/20)*std(im(:));
y = y + (randn(size(y)) + 1i*randn(size(y)))*sigma_noise/sqrt(2);
% Display the downsampled image
imagesc_gray(real(At(y)), 2 , 'Measured image',121);
% Tolerance on noise
epsilon = sqrt(chi2inv(0.99, 2*numel(ind))/2)*sigma_noise;
% Solve TVDN problem
sol = solve_tvdn(y, epsilon, A, At, param);
% Show reconstructed image
%%
imagesc_gray(abs(sol), 2 , 'Reconstructed image',122);

%%
close_unlocbox();

