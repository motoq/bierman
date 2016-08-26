function [Q, R] = mth_qr(A)
% MTH_QR decomposes an input covariance A matrix using "thin" QR decomposition
% using the Gram-Schmidt method as described in "Optimal Estimation of Dynamic
% Systems" by Crassidis and Junkins instead of the build in Octave/Matlab
% version that uses Householder transformations (this is just for comparison
% purposes).
% A = Q*R
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
%   A   [MxN] Input matrix
% Return:
%   Q   [MxN] orthonormal matrix
%   R   [NxN] upper triangular matrix
%
% Author:  Kurt Motekew    20160816
% 

  [~, n] = size(A);
  R = zeros(n);
  Q = A;
  for kk = 1:n
    qk = Q(:,kk);
    rkk = norm(qk);
    R(kk,kk) = rkk;
    qk = qk/rkk;
    Q(:,kk) = qk;
    for jj = (kk+1):n
      qk = Q(:,kk);
      qj = Q(:,jj);
      R(kk,jj) = dot(qk,qj);
      qk = qk*R(kk,jj);
      qj = qj - qk;
      Q(:,jj) = qj;
    end
  end
