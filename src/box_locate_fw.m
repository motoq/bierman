function [phat, SigmaP, itr] = box_locate_fw(tkr_pos, y, W)
% BOX_LOCATE_FW Geolocates a tracked object within a boxed volume given
% range only tracker locations, measurements, and a measurement weighting
% matrix using a WLS solution.  Observation errors are not assumed to be
% uncorrelated, leading to use of a full weighting matrix without the
% measurement accumulation allowed with a block diagonal structure.
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
%   tkr_pos   A [3xM] matrix of M tracker locations
%   y         Array of M distance measurements
%   W         Range uncertianty weighting matrix, [MxM]
%
% Return:
%   phat     Estimated location, [3x1]
%   SigmaP   Location covariance, [3x3]
%   itr      Number of iterations
%
% Kurt Motekew   2016/08/13
%
  maxitr = 50;
  tol = .0000001;

  nmeas = size(tkr_pos, 2);

  phat = [0.5 0.5 0.5]';
  r = zeros(nmeas,1);
  Ap = zeros(nmeas, 3);
    % Populate residual and partial derivative matrices
  for itr = 1:maxitr
    for ii = 1:nmeas
      sc = phat - tkr_pos(:,ii);
      yc = norm(sc);
      r(ii) = y(ii) - yc;
      Ap(ii,:) = est_drng_dloc(tkr_pos(:,ii), phat);
    end
    ApTWAp = Ap'*W*Ap;
    SigmaP = ApTWAp^-1;
    dp = SigmaP*Ap'*W*r;
    phat = phat + dp;
    if norm(dp) < tol
      break;
    end
  end

