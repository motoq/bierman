function L2 = mth_chol_upd_l(L, beta, v)
% MTH_CHOL_UPD Given L = sqrt(P) such that P = L*L', compute L2 such
% that P2 = P + c*a*a' where 'c' is a scalar, 'a' is a vector, and
% P2 = L2*L2'.
%
%-----------------------------------------------------------------------
% Copyright 2019 Kurt Motekew
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%-----------------------------------------------------------------------
%
% Inputs:
%   L   Lower triangular square root of an [NxN] matrix P such that
%       P = L*L'
%   c   Scalar
%   a   Vector [Nx1]
% Return:
%   L   Updated lower triangular square root such that
%       L = sqrt(L*L' +  c*a*a')
%
% Author:  Kurt Motekew    20190110
%
% Ref:  Krause & Igel, "A More Efficient Rank-one Covariance Matrix Update for
%       Evolution Strategies"
% 

  w = v;
  b = 1;
  n = size(L,1);
  L2 = zeros(n);
  for jj = 1:n
    L2(jj,jj) = sqrt(L(jj,jj)*L(jj,jj) + beta*w(jj)*w(jj)/b);
    gamma = L(jj,jj)*L(jj,jj)*b + beta*w(jj)*w(jj);
    for kk = (jj+1):n
      w(kk) =  w(kk) - w(jj)*L(kk,jj)/L(jj,jj);
      L2(kk,jj) = L2(jj,jj)*L(kk,jj)/L(jj,jj) +...
                  L2(jj,jj)*beta*w(jj)*w(kk)/gamma;
    end
    b = b + beta*w(jj)*w(jj)/(L(jj,jj)*L(jj,jj));
  end
