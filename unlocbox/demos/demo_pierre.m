%DEMO_PIERRE Demo to solve a particular l1 l2 problem
%
%   The problem can be expressed like this
%
%        argmin_{c,b} ||s - Psi c - Phi b||_2^2 + mu1 ||c||_1 + mu2 ||b||_1
%
%   Where s are the measurements, Psi the Fourier matrix and
%   Phi=Phi*M with M a diagonal matrix with +1,-1 random values.
%   
%   We will use generalized forward backward to solve this problem. The
%   gradients of 
%
%      `||s - Psi c - Phi b||_2^2 
%
%   are
%
%      nabla_{c}f(c,b) = 2 Psi^* (Psi c + Phi b - s)
%
%      nabla_{b}f(c,b) = 2 Phi^* (Psi c + Phi b - s)
%
%   In this code the variable b and c will be stack into one single
%   vector of size 2N
%
%   Results
%   -------
%
%   Figure 1: Results of the reconstruction
%
%      The support of the signal is recovered.
%
%   Url: http://lts2research.epfl.ch/unlocbox/doc/demos/demo_pierre.php

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
% Date: 2 October 2013

%% Initialisation

clear all;
close all;

% Loading toolbox
init_unlocbox();

verbose = 1;    % verbosity level

%% General parameter for the problem
mu1=0.001;
mu2=0.001;

%% Creation of the Problem

% Size of the problem;
N=1024;
K=50;

% Creation of the matrix M
M=rand(N,1);
M=M>0.5;
M=M*2-1;
% We keep M in a vector form.
%M= diag(M);

% Creation of the operator Psi
Psi=@(x) 1/sqrt(N)*fft(x);
Psit=@(x) sqrt(N)*ifft(x);
nu_Psi=1; % Frame bound

% Creation of th operator Phi
Phi= @(x) M.*Psi(x);
Phit= @(x) Psit(M.*x);
nu_Phi=1; % Frame bound


% Create a K sparse signal
sigb = zeros(N, 1); I = randperm(N);
sigb(I(1:K)) = randn(K, 1);
sigb = sigb/norm(sigb);

% Create a K sparse signal
sigc = zeros(N, 1); I = randperm(N);
sigc(I(1:K)) = randn(K, 1);
sigc = sigc/norm(sigc);


% Creation of the measurements
s=Psi(sigc)+Phi(sigb);

%% Define the UNLocBoX proximity operators and gradient

%Define mask functions
Mc = @(x) x(1:N);         % select the variable c
Mb = @(x) x(N+1:end);     % select the variable b
Cc = @(x,c) [c;Mb(x)];    % recreate all opt variable from c
Cb = @(x,b) [Mc(x);b];    % recreate all opt variable from b

% Functional mu1*||c||_1
param_l1_1.verbose = verbose - 1;
f1.prox = @(x,T) Cc(x,prox_l1(Mc(x),T*mu1,param_l1_1));
f1.eval = @(x)   mu1*norm(Mc(x),1);

% Functional mu2*||c||_1
param_l1_2.verbose = verbose - 1;
f2.prox = @(x,T) Cb(x,prox_l1(Mb(x),T*mu1,param_l1_2));
f2.eval = @(x)   mu2*norm(Mb(x),1);

% Define the gradient of ||s - Psi c - Phi b||_2^2 
g.grad = @(x)  [ Psit( Psi(Mc(x)) + Phi(Mb(x)) - s ) ; ...
                 Phit( Psi(Mc(x)) + Phi(Mb(x)) - s ) ];
g.eval = @(x)    norm(s-(Psi(Mc(x))+Phi(Mb(x))));


%% Define solver parameter

paramfb.maxit=1000;
paramfb.gamma=1/(nu_Phi+nu_Psi);
paramfb.verbose=verbose;

% Solve the problem
x0= eps*ones(2*N,1);
sol=generalized_forward_backward(x0,{f1,f2},g,paramfb);


%% Display the results
figure(1);
subplot(211)
t=1:N;
plot(t,abs(sigb),'xb',t,abs(Mb(sol)),'or');
legend('Original b','Recovered b');
subplot(212)
plot(t,abs(sigc),'xb',t,abs(Mc(sol)),'or');
legend('Original c','Recovered c');

