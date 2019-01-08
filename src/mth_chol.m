function R = mth_chol(P)
% MTH_CHOL decomposes an input covariance P into an upper triangular matrix
% R via Cholesky decomposition such that P = R'*R.  Note the difference
% between mth_sqrtm.
%
%-----------------------------------------------------------------------
% Copyright 2018 Kurt Motekew
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%-----------------------------------------------------------------------
%
% Inputs:
%   P   Positive definite symmetrix covariance, [NxN]
% Return:
%   R   Square root of P, upper triangular, P = R'*R
%
% Author:  Kurt Motekew    20181207
% 
% Ref:  Trefethen & Bau, Numerical Linear Algebra
%       SIAM
%

  m = size(P,1);
  R = P;
  for kk = 1:m
    for jj = (kk+1):m
      R(jj,jj:m) = R(jj,jj:m) - R(kk,jj:m)*R(kk,jj)/R(kk,kk);
    end
    R(kk,kk:m) = R(kk,kk:m)/sqrt(R(kk,kk));
  end

  for kk = 2:m
    jj = 1:(kk-1)
      R(kk,jj) = 0;
    end
  end

