function [U, D] = mth_udut2(P)
% MTH_UDUT2 decomposes an input covariance P into U and D such that U is a unit
% upper diagonal matrix and D is a diaginal matrix with P = U*D*U'.
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
%   U   Unit upper diagonal triangular matrix [NxN]
%   D   Diagonal matrix [NxN]
%
% Note for conversion to other languages:  P is reused internally.  Since
% Matlab and Octave pass by copy, it doesn't matter if P is destroyed.
%
% Author:  Kurt Motekew    20160802
% 
  n = size(P,2);
  U = eye(n);
  D = zeros(n);
  for jj = n:-1:2
    D(jj,jj) = P(jj,jj);
    alpha = 1/D(jj,jj);
    for kk = 1:(jj-1)
      beta = P(kk,jj);
      U(kk,jj) = alpha*beta;
      for ii = 1:kk
        P(ii,kk) = P(ii,kk) - beta*U(ii,jj);
      end
    end
  end
  D(1,1) = P(1,1);

