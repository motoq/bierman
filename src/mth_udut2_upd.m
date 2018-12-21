function [U2, D2] = mth_udut2_upd(U, D, c, a)
% MTH_UDUT2_UPD Given P = UDU' where U is a unit diagonal triangular
% matrix and D is a diagonal matrix, computes a new U and D for
% P2 = P + c*a*a' where 'c' is a scalar and 'a' is a vector.
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
%   U   Unit upper diagonal triangular matrix [NxN]
%   D   Diagonal matrix [NxN]
%   c   Scalar
%   a   Vector [Nx1]
% Return:
%   U2   Updated unit upper diagonal triangular matrix [NxN]
%   D2   Updated diagonal matrix [NxN]
%
% Author:  Kurt Motekew    20181127
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

