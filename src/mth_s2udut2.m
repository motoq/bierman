function [U, D] = mth_s2udut2(S)
% MTH_S2UDUT2 Converts an upper triangular square root matrix into
% a unit upper triangular matrix and a diagonal matrix such that
% U*D*U' = S*S'
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
%   S   Upper triangular square root of an [NxN] matrix P such that
%       P = S*S' and S = U*sqrt(D)
% Return:
%   U   Unit diagonal upper triangular matrix [NxN]
%   D   Diagonal matrix [NxN]
%
% Author:  Kurt Motekew    20181220
% 

  n = size(S,1);
  U = eye(n);
  Dsr = diag(S);                                 % Sqrtm(D)

  for ii = (n-1):-1:1
    for jj = n:-1:(ii+1)
      U(ii,jj) = S(ii,jj)/Dsr(jj);
    end
  end
  D = diag(Dsr.*Dsr);
