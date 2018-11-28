function [U2, D2] = mth_udut2_upd(U, D, c, a)
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
  n = size(D,2);
  U2 = eye(n);
  D2 = zeros(n);

  for jj = n:-1:2
    s = a(jj);
    d = D(jj,jj) + c*s*s;
    D2(jj,jj) = d;
    b = c/d;
    beta = s*b;
    c = b*D(jj,jj);
    for ii = (jj-1):-1:1
      a(ii) = a(ii) - s*U(ii,jj);
      U2(ii,jj) = U(ii,jj) + beta*a(ii);
    end
  end
  D2(1,1) = D(1,1) + c*a(1)*a(1);

