function ai = mth_lpoly_fit(x, y, order)
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
%   x  [1xM]
%   y  [1xM]
%
% Kurt Motekew   2018/11/04
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
  for ii = 1:m
    Ai(1) = 1;
    Ai(2) = x(ii);
    for jj = 3:n
      Ai(jj) = x(ii)*Ai(jj-1);
    end
    ATA = ATA + Ai'*Ai;
    ATy = ATy + Ai'*y(ii);
  end

  ai = ATA^-1*ATy;
  ai = ai';

