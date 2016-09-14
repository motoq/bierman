function [x, S] = est_upd_potter(x, S, A, delta, sw)
% EST_UPD_POTTER Updates the apriori estimate and covariance via the Potter
% mechanization of the Kalman filter given a single new observation.
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
%   S         A priori estimate covariance square root (P = SS'), [NxN]
%   A         Partial of obs w.r.t. x, [1xN]
%   delta     Observation residual, scalar
%   sw        Inverse of observation uncertainty, 1/sigma, scalar
%
% Return:
%   x        Updated estimated, [Nx1]
%   S        Updated covariance square root (P = SS'), [NxN]
%
% Kurt Motekew   2016/08/04
%
%
% Ref:  G. J. Bierman, Factorization Methods for
%       Discrete Sequential Estimation, Dover Publications, Inc.,
%       Mineola, NY, 1977, pp. 30-31.
%

    % Normalize, Covariance obs = I
  delta = sw*delta;                              % Scaled predicted residual
  A = sw*A;                                      % Scaled partials
  
  vtrans = A*S;                                  % [1xN]
  sigma = 1/(vtrans*vtrans' +1);                 % Predicted residual covariance
  K = S*vtrans';                                 % Kalman gain
  x = x + K*(delta*sigma);                       % State update
  lambda = sigma/(1 + sqrt(sigma));
  S = S - (lambda*K)*vtrans;                     % Stabilized covariance update
