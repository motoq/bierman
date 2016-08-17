function [phat, SigmaP, R, z, itr] = box_locate_hh(tkr_pos, y, SqrtW)
% BOX_LOCATE_HH Geolocates a tracked object within a boxed volume given
% range only tracker locations, measurements, and a measurement weighting
% matrix using Householder triangularization.
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
%   SqrtW     Square root of range uncertianty weighting matrix, [MxM]
%
% Return:
%   phat     Estimated location, [3x1]
%   SigmaP   Location covariance, [3x3]
%   R
%   z
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
  for itr = 1:maxitr
      % Build residual and partials matrices
    for ii = 1:nmeas
      sc = phat - tkr_pos(:,ii);
      yc = norm(sc);
      r(ii) = y(ii) - yc;
      Ap(ii,:) = est_drng_dloc(tkr_pos(:,ii), phat);
    end
      % Form information array and decompose for solution
    InfoArray = [SqrtW*Ap SqrtW*r];
    Rz = mth_householder_tri(InfoArray, 1);
    R = Rz(1:3,1:3);
    z = Rz(1:3,4);
    dp = mth_trisol(R, z);
    phat = phat + dp;
    if norm(dp) < tol
      break;
    end
  end
  Rinv = mth_triinv(R);
  SigmaP = Rinv*Rinv';

