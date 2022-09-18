function y = tv_normnd(u,weights)
%TV_NORMND N Dimentional TV norm
%   Usage:  tv_normnd(x,weights)
%
%   Input parameters:
%         x     : Input data (N dimentional matrix)
%         weights: Weights
%   Output parameters:
%         sol   : Norm
%
%   Compute the N-dimentional TV norm of x
%
%   Url: http://lts2research.epfl.ch/unlocbox/doc/misc/tv_normnd.php

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


warning('This function will be deleted in future release of matlab, use norm_ndtv instead')

sz = size(u);
        anisotropic = 1;
        dim = length(sz);
        temp = zeros(sz);
        if nargin<2
            weights = ones(dim,1);
        end
    if anisotropic == 1
        for d = 1:dim
            sz_temp = sz;
            sz_temp(d) = 1;
            tv(d).grad = weights(d)*cat(d,diff(u,1,d),zeros(sz_temp));
            temp = abs(tv(d).grad)+ temp;
        end
        y = sum(temp(:));
    else
        for d = 1:dim
            sz_temp = sz;
            sz_temp(d) = 1;
            tv(d).grad = weights(d)*cat(d,diff(u,1,d),zeros(sz_temp));
            temp = abs(tv(d).grad.^2)+ temp;
        end
        temp = sqrt(temp);
        y = sum(temp(:));
    end
end
