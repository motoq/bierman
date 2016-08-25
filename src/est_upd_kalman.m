function [phat, SigmaP] = est_upd_kalman(phat0, P0, Ap, r, sw)
% EST_UPD_KALMAN Updates the apriori estimate and covariance via a stabilized
% Kalman filter given a single new observation.
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
%   phat0     A priori estimate, [Nx1]
%   P0        A priori estimate covariance, [NxN]
%   Ap        Partial of obs w.r.t. solve phat0, [1xN]
%   r         Observation residual, scalar
%   sw        Inverse of observation uncertainty, 1/sigma, scalar
%
% Return:
%   phat     Updated estimated, [Nx1]
%   SigmaP   Updated covariance, [NxN]
%
% Kurt Motekew   2016/08/04
%
%
% Ref:  G. J. Bierman, Factorization Methods for
%       Discrete Sequential Estimation, Dover Publications, Inc.,
%       Mineola, NY, 1977, pp. 28-30.
%

    % Normalize, Covariance obs = I
  delta = sw*r;                                  % Scaled predicted residual
  Ap = sw*Ap;                                    % Scaled partials
  
  v = P0*Ap';
  sigma = Ap*v +1;                               % Predicted residual covariance
  K = v/sigma;                                   % Kalman gain
  phat = phat0 + K*delta;                        % State update
  Pbar = P0 - K*v';                              % Optimal covariance update
  v = Pbar*Ap';
  SigmaP = (Pbar - v*K') + K*K';                 % Stabilized covariance update
