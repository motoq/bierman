function x = mth_trisol(U, y)
% MTH_TRISOL Uses backwards substitution to solve x = U^-1 * y given U
% is upper triangular.
%
%-----------------------------------------------------------------------
% Copyright 2016 Kurt Motekew
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%-----------------------------------------------------------------------
%
% Inputs:
%   U   Upper triangular [MxM] matrix.
%   y   [Mx1] vector
%
% Return:
%   x   U^-1 * y
%
% Author:  Kurt Motekew    20160813
% 

  n = size(U,1);
  x = zeros(n,1);
  for jj = n:-1:1
    tmp = 0;
    for kk = (jj+1):n
      tmp = tmp + U(jj,kk)*x(kk);
    end
    x(jj) = (y(jj) - tmp)/U(jj,jj);
  end
