function x = mth_lpoly_eval(ai, t)
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
%   ai  [1xM]
%   t   scalar
%
% Kurt Motekew   2018/11/05
%

    % Reduce order to maximum supported by observations
  n = size(ai,2);
  m = size(t, 2);

  x = zeros(1,m);
  x(1,:) = ai(1);
  ti = ones(1,m);
  for ii = 2:n
    ti = ti.*t;
    x = x + ai(ii).*ti;
  end
