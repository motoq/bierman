function Uinv = mth_triinv(U)
% MTH_TRIINV Inverts an upper triangular matrix
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
%
% Return:
%   Uinv   Invers of input U
%
% Author:  Kurt Motekew    20160813
% 

  n = size(U,1);
  Uinv = zeros(n);

  Uinv(1,1) = 1/U(1,1);
  for jj = 2:n
    Uinv(jj,jj) = 1/U(jj,jj);
    jm1 = jj - 1;
    for kk = 1:jm1
      tmp = 0;
      for ii = kk:jm1
        tmp = tmp - Uinv(kk,ii)*U(ii,jj);
      end
      Uinv(kk,jj) = tmp*Uinv(jj,jj);
    end
  end
