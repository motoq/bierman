function [x, Px, Pxy] = est_upd_schmidt(x, Px, Ax, Py, Ay, Pxy, delta, Pz)
% EST_UPD_SCHMIDT Updates the apriori estimate and covariance via the Schmidt-
% Kalman "minimum covariance" Kalman filter.
% mechanization of the Kalman filter.
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
%   x         A priori estimate, [Nx1]
%   Px        A priori estimate covariance, [NxN]
%   Ax        Partial of obs w.r.t. x, [MxN]
%   Py        Consider parameter covariance, [NyxNy]
%   Ay        Partial of obs w.r.t. consider params, [MxNy]
%   Pxy       Estimate/consider corss covariance, [NxNy]
%   delta     Observation residual, scalar
%   Pz        Observation covariance, [MxM]
%
% Return:
%   x        Updated estimated, [Nx1]
%   Px       Updated estimate covariance, [NxN]
%   Pxy      Updated estimate/consider cross covariance, [NxNy]
%
% Kurt Motekew   2017/02/14
%
%
% Ref:  Woodbury, D.P., Junkins, J.L.,
%       On the Consider Kalman Filter
%       AIAA 2010-7752
%

  K = (Px*Ax' + Pxy*Ay')*...
      (Ax*Px*Ax' + Ax*Pxy*Ay' + Ay*Pxy'*Ax' + Ay*Py*Ay' + Pz)^-1;
  x = x + K*delta;
  n = size(x,1);
  JF = eye(n) - K*Ax;
  Px  = JF*Px  - K*Ay*Pxy';
  Pxy = JF*Pxy - K*Ay*Py;

