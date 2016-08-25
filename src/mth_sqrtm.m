function [U] = mth_sqrtm(P)
% MTH_SQRTM decomposes an input covariance P into an upper triangular matrix
% U via Cholesky decomposition such that P = U*U'.
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
%   P   Positive definite symmetrix covariance, [NxN]
% Return:
%   U   Square root of P, upper triangular, P = U*U'
%
% Author:  Kurt Motekew    20160824
% 
% Ref:  G. J. Bierman, Factorization Methods for
%       Discrete Sequential Estimation, Dover Publications, Inc.,
%       Mineola, NY, 1977, p. 53.
%

  n = size(P,2);
  U = zeros(n);

  for jj = n:-1:2
    U(jj,jj) = sqrt(P(jj,jj));
    alpha = 1/U(jj,jj);
    for kk = 1:(jj-1)
      U(kk,jj) = alpha*P(kk,jj);
      beta = U(kk,jj);
      for ii = 1:kk
        P(ii,kk) = P(ii,kk) - beta*U(ii,jj);
      end
    end
  end
  U(1,1) = sqrt(P(1,1));

