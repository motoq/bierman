function [phat, SigmaP, itr] = box_locate(tkr_pos, y, W)
% BOX_LOCATE Geolocates a tracked object within a boxed volume given
% range only tracker locations, measurements, and a measurement weighting
% matrix.
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
%   tkr_pos   A [3xM] matrix of N tracker locations
%   y         Array of M distance measurements
%   W         Range uncertianty weighting matrix, 1/sigma_range^2, scalar
%
% Return:
%   phat     Estimated location, [3x1]
%   SigmaP   Location covariance, [3x3]
%   itr      Number of iterations
%
% Kurt Motekew   2014/10/25
%
  maxitr = 50;
  tol = .0000001;

  nmeas = size(tkr_pos, 2);

  phat = [0.5 0.5 0.5]';
  yc(nmeas) = 0;
  for itr = 1:maxitr
    for ii = 1:nmeas
      sc = phat - tkr_pos(:,ii);
      yc(ii) = norm(sc);
    end
    r = y - yc;
    ApTWAp = zeros(3);
    ApTWr = zeros(3,1);
      % Accumulate measurements
    for ii = 1:nmeas
      Api = est_drng_dloc(tkr_pos(:,ii), phat);
      ApTWAp = ApTWAp + Api'*W*Api;
      ApTWr = ApTWr + Api'*W*r(ii);
    end
    SigmaP = ApTWAp^-1;
    dp = SigmaP*ApTWr;
    phat = phat + dp;
    if norm(dp) < tol
      break;
    end
  end

