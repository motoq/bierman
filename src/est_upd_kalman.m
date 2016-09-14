function [x, P] = est_upd_kalman(x, P, A, delta, sw)
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
%   x         A priori estimate, [Nx1]
%   P         A priori estimate covariance, [NxN]
%   A         Partial of obs w.r.t. solve x, [1xN]
%   delta     Observation residual, scalar
%   sw        Inverse of observation uncertainty, 1/sigma, scalar
%
% Return:
%   x        Updated estimated, [Nx1]
%   P        Updated covariance, [NxN]
%
% Kurt Motekew   2016/08/04
%
%
% Ref:  G. J. Bierman, Factorization Methods for
%       Discrete Sequential Estimation, Dover Publications, Inc.,
%       Mineola, NY, 1977, pp. 28-30.
%

    % Normalize, Covariance obs = I
  delta = sw*delta;                              % Scaled predicted residual
  A = sw*A;                                      % Scaled partials
  
  v = P*A';
  sigma = A*v + 1;                               % Predicted residual covariance
  K = v/sigma;                                   % Kalman gain
  x = x + K*delta;                               % State update
  P = P - K*v';                                  % Optimal covariance update
  v = P*A';
  P = (P - v*K') + K*K';                         % Stabilized covariance update
