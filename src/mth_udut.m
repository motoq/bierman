function [UD] = mth_udut(P)
% MTH_UDUT decomposes an input covariance P into U and D such that U is a unit
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
%   UD  Combined unit diagonal upper triangular matrix U and diagonal matrix
%       D.  D is along the diagonal since U's diagonal would be 1's.  [NxN]
%
% Note for conversion to other languages:  P is reused internally.  Since
% Matlab and Octave pass by copy, it doesn't matter if P is destroyed.
%
% Author:  Kurt Motekew    20160802
% 
  n = size(P,2);
  UD = zeros(n);
  D = zeros(1,n);
  for jj = n:-1:2
    D(jj) = P(jj,jj);
    alpha = 1/D(jj);
    for kk = 1:(jj-1)
      beta = P(kk,jj);
      UD(kk,jj) = alpha*beta;
      for ii = 1:kk
        P(ii,kk) = P(ii,kk) - beta*UD(ii,jj);
      end
    end
  end
  D(1) = P(1,1);

  for ii = 1:n
    UD(ii,ii) = D(ii);
  end

  

