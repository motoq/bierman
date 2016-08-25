function [phat, S] = est_upd_potter(phat0, S0, Ap, r, sw)
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
%   phat0     A priori estimate, [Nx1]
%   S0        A priori estimate covariance square root (P = SS'), [NxN]
%   Ap        Partial of obs w.r.t. phat0, [1xN]
%   r         Observation residual, scalar
%   sw        Inverse of observation uncertainty, 1/sigma, scalar
%
% Return:
%   phat     Updated estimated, [Nx1]
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
  delta = sw*r;                                  % Scaled predicted residual
  Ap = sw*Ap;                                    % Scaled partials
  
  vtrans = Ap*S0;                                % [1xN]
  sigma = 1/(vtrans*vtrans' +1);                 % Predicted residual covariance
  K = S0*vtrans';                                % Kalman gain
  phat = phat0 + K*(delta*sigma);                % State update
  lambda = sigma/(1 + sqrt(sigma));
  S = S0 - (lambda*K)*vtrans;                    % Stabilized covariance update
