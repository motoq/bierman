function [x_hat, S_hat] = est_upd_srukf(x_bar, S_bar, Chi, w_m, sr_w_c,...
                                        Y, y, Sr_Rn)
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
%   S_bar  Estimate covariance, square root, currently LOWER triangular [mXm]
%   Chi    Sigma vectors (used to form x_bar and P_bar), [mXn]
%          where n is the number of sigma vectors.
%   w_m    Estimate weighting factors, [1Xn]
%   w_c    Covariance weighting factors, [1Xn]
%   Y      Sigma vector based computed observations, [num_obs X n]
%   y      Observations, [num_obx X 1]
%   Rn     Process noise, [mxm]
%
% Return:
%   x_hat  State estimate update based on observations
%   S_hat  Updated estimate covariance square root, upper triangular, [mxm]
%
% Kurt Motekew   2018/11/14
%
%

  P_bar = S_bar*S_bar';
  w_c = sr_w_c.*sr_w_c;

  dim = size(Chi,1);
  n_sigma_vec = size(Chi, 2);
  n_obs = size(Y,1);

    % Computed obs based on propagated state
  y_bar = zeros(n_obs,1);
  for kk = 1:n_sigma_vec
    y_bar = y_bar + w_m(kk)*Y(:,kk);
  end
    % Observation update
  SigmaXY = zeros(dim,n_obs);
  AT = [Y  Sr_Rn];
  for kk = 1:n_sigma_vec     
    y_minus_ybar = Y(:,kk) - y_bar;
    chi_minus_xbar = Chi(:,kk) - x_bar;
    AT(:,kk) = sr_w_c(kk)*y_minus_ybar;
    SigmaXY = SigmaXY + w_c(kk)*(chi_minus_xbar*y_minus_ybar');
  end
  [~, S_Y_bar] = mth_qr(AT');
  SigmaY_bar = S_Y_bar'*S_Y_bar;

  K = SigmaXY*SigmaY_bar^-1;
  x_hat = x_bar + K*(y - y_bar);

  U = K*S_Y_bar';
  n = size(U,2);

  S_hat = mth_sqrtm(P_bar);
  for kk = 1:n
    S_hat = mth_chol_upd(S_hat, -1, U(:,kk));
  end

