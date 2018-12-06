function x = mth_lpoly_eval(ai, t)
% MTH_LPOLY_EVAL Evaluates a Lagrange polynomial.
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
%   ai  Polynomial coefficients, [1xN+1].  An nth order polynomial will
%       contain n + 1 coefficients:
%         y = ai(1) + ai(2)*x } ai(3)x^2 + ... + ai(n+1)*x^n
%   t   Independent parameter, [1XM]
%
% Return:
%   x   Function value, [1XM]
%
% Kurt Motekew   2018/12/05
%

    % n = order + 1
  n = size(ai,2);
    % m = number of function points to evaluate
  m = size(t, 2);

  x = zeros(1,m);
  x(1,:) = ai(1);
  ti = ones(1,m);
  for ii = 2:n
    ti = ti.*t;
    x = x + ai(ii).*ti;
  end
