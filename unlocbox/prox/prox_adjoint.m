function sol  = prox_adjoint( x,gamma,f )
%PROX_ADJOINT Proximal operator of the adjoint function of f
%   Usage:   sol=prox_adjoint(x, gamma, f);
%
%   Input parameters:
%         x     : Input signal.
%         gamma : Regularization parameter.
%         f     : Function
%   Output parameters:
%         sol   : Solution.
%         infos : Structure summarizing informations at convergence
%
%    PROX_ADJOINT( x,gamma,f ) solves:
%
%      sol = argmin_{z} 0.5*||x - z||_2^2 + gamma * f*
%
%   where f^ is the adjoint of f. This problem is solved thanks to the
%   Moreau's identity.
%
%   Warning: f needs to be a proper convex lower semi continuous function.
%
%
%   Url: http://lts2research.epfl.ch/unlocbox/doc/prox/prox_adjoint.php

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
% Date: 31 May 2013
%

if gamma==0
    warning(['gamma = 0. This is problably not correct.' ...
             ' We replace it by eps to keep going.']);
    gamma=eps;
end

% Optional input arguments
if nargin<3,     
    error('Two few input arguments!');
end

sol_a = f.prox(x/gamma,1/gamma);

sol= x-sol_a;


end


