function [x_hat, P_hat] = est_upd_ukf(x_bar, P_bar, Chi, w_m, w_c,...
                                      Y, y, Rn)
% EST_UPD_UKF Given the current estimate and covariance, update with
% sigma vector based parameters and computed observations.
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
%   x_bar  Current estimate [mX1]
%   P_bar  Estimate covariance [mXm]
%   Chi    Sigma vectors (used to form x_bar and P_bar), [mXn]
%          where n is the number of sigma vectors.
%   w_m    Estimate weighting factors, [1Xn]
%   w_c    Covariance weighting factors, [1Xn]
%   Y      Sigma vector based computed observations, [num_obs X n]
%   y      Observations, [num_obs X 1]
%   Rn     Observation covariance, [num_obs X num_obs]
%
% Return:
%   x_hat  State estimate update based on observations
%   P_hat  Updated estimate covariance
%
% Kurt Motekew   2018/11/14
%
% Ref:  Rudolph van der Merwe & Eric A. Wan, "The Unscented Kalman Filter
%       for Nonlinear Estimation"
%

  dim = size(Chi,1);
  n_sigma_vec = size(Chi, 2);
  n_obs = size(Y,1);

    % Computed obs based on propagated state
  y_bar = zeros(n_obs,1);
  for kk = 1:n_sigma_vec
    y_bar = y_bar + w_m(kk)*Y(:,kk);
  end
    % Observation update
  SigmaY_bar = zeros(n_obs);
  SigmaXY = zeros(dim,n_obs);
  for kk = 1:n_sigma_vec     
    y_minus_ybar = Y(:,kk) - y_bar;
    chi_minus_xbar = Chi(:,kk) - x_bar;
    SigmaY_bar = SigmaY_bar + w_c(kk)*(y_minus_ybar*y_minus_ybar');
    SigmaXY = SigmaXY + w_c(kk)*(chi_minus_xbar*y_minus_ybar');
  end
  SigmaY_bar = SigmaY_bar + Rn;
  K = SigmaXY/SigmaY_bar;
  x_hat = x_bar + K*(y - y_bar);
  P_hat = P_bar - K*SigmaY_bar*K';
