function Q = mth_hh_q(A)
% MTH_HH_Q Generates the Q portion of the full QR decomposition of a
% matrix via Householder transformations.
%
%-----------------------------------------------------------------------
% Copyright 2022 Kurt Motekew
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%-----------------------------------------------------------------------
%
% Inputs
%   A  MxN matrix, M>=N for a Q that offers full triangularization
%
% Return
%  Q  Orthonormal matrix such that Q'A is upper triangular
%
% Kurt Motekew  2022/11/10
%

  [mm,nn] = size(A);
    % Starting transformation to build on
  Q = eye(mm);

    % For each iteration, the upper-left portion of each matrix/vector
    % is untouched (allowing for full QR to transform in place).
    % Separate x, v, and u vectors are used for clarity vs.
    % speed/memory (wouldn't be using Matlab in the first place if
    % speed/memory were of concern...).
  T = zeros(mm,mm);                         % Householder reflection
  x = zeros(mm,1);                          % Current matrix column
  v = zeros(mm,1);                          % Current reflection normal 
  u = zeros(mm,1);                          % Current unit reflection normal
  e1 = zeros(mm,1);                         % Shrinking ihat basis vector

    % Q = T1*T2...*Tn with the upper left of T becoming an increasingly
    % large identity matrix and T becoming increasingly (block) diagonal.
  for ii = 1:nn
    e1(ii) = 1;
    x(ii:mm) = A(ii:mm,ii);
    v(ii:mm) = x(ii:mm) + sign(x(ii))*norm(x(ii:mm))*e1(ii:mm);
    u(ii:mm) = v(ii:mm)/norm(v(ii:mm));
    T(ii:mm,ii:mm) = eye(mm-ii+1) - 2*u(ii:mm)*u(ii:mm)';
    A(ii:mm,ii:nn) = T(ii:mm,ii:mm)*A(ii:mm,ii:nn); 
    Q = Q*[eye(ii-1) zeros(ii-1,mm-ii+1) ; zeros(mm-ii+1,ii-1) T(ii:mm,ii:mm)];
  end
