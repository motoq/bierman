function K = mth_tri_inv_mat_mul(P, S)
% MTH_TRI_INV_MAT_MUL Uses forward and backwards substitution to solve
% for K given P and S such that K*S'*S = P and S is upper triangular.
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
%   P   [MxN] Matrix
%   S   [NxN] Upper triangular matrix
%
% Return:
%   K   P*S^-1*(S^1)', [MxN]
%
% Author:  Kurt Motekew    20190122
% 

  [m, n] = size(P);

    % Solve for K*S'
  Kst = zeros(m,n);
  for ii = 1:m
    Kst(ii,1) = P(ii,1)/S(1,1);
    for jj = 2:n
      s = P(ii,jj);
      for kk = jj:-1:2
        s = s - Kst(ii,kk-1)*S(kk-1,jj);
      end
      Kst(ii,jj) = s/S(jj,jj);
    end
  end

    % Now find K
  K = zeros(m,n);
  for ii = 1:m
    K(ii,n) = Kst(ii,n)/S(n,n);
    for jj = (n-1):-1:1
      s = Kst(ii,jj);
      for kk = n:-1:(jj+1)
        s = s - K(ii,kk)*S(jj,kk);
      end
      K(ii,jj) = s/S(jj,jj);
    end
  end

