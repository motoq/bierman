function ai = mth_lpoly_fit(x, y, order)
% MTH_LPOLY_FIT Creates Lagrange polynomial coefficients given coordinate
% points.  The system can be overdetermined, resulting in a least squares
% fit.
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
%   x      Independent variable values [1xM]
%   y      Observed values [1xM]
%   order  Polynomial order N, N <= M
%
% Output
%   Polynomial coefficients, [1xN+1].  An nth order polynomial will
%   contain n + 1 coefficients:
%     y = ai(1) + ai(2)*x } ai(3)x^2 + ... + ai(n+1)*x^n
%
% Kurt Motekew   2018/12/05
%

    % Reduce order to maximum supported by observations
  m = size(x,2);
  if nargin == 3
    if m >= (order + 1)
      n = order + 1;
    else
      n = m;
    end
  else
    n = m;
  end
   
    % nObs X nCoeff
  Ai = zeros(1,n);
  ATA = zeros(n);
  ATy = zeros(n,1);
    % Accumulate obs
  for ii = 1:m
    Ai(1) = 1;
    Ai(2) = x(ii);
    for jj = 3:n
      Ai(jj) = x(ii)*Ai(jj-1);
    end
    ATA = ATA + Ai'*Ai;
    ATy = ATy + Ai'*y(ii);
  end

    % Solve for coefficients via QR decomp and backwards substitution
  [Q, R] = mth_qr(ATA);
  qtaty = Q'*ATy;
  ai = mth_trisol(R, qtaty);
    % Return as row vector
  ai = ai';

