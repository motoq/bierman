function [x_bar, S_bar] = est_pred_srukf(Chi, w_m, sr_w_c, Sr_Rv)
% EST_PRED_UKF Given propagated sigma vectors, updates an estimate's
% predicted state and covariance.
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
%   Chi      Sigma vectors to be condensed into a new estimate, [mxn]
%            where m is the dimension of the state vector and n is the
%            number of sigma vectors.
%   w_m      Estimate weighting factors, [1xn]
%   sr_w_c   Square root of covariance weighting factors, [1xn]
%   Sr_Rv    Square root of process noise, [mxm]
%
% Return:
%   x_bar  Updated estimate based on propagated sigma vectors, [mx1]
%   S_bar  Updated estimate covariance square root based on propagated
%          sigma vectors, LOWER triangular for now [mxm]
%
% Kurt Motekew   2018/12/19  Base4d on  est_pred_ukf
%
%

  dim = size(Chi,1);
  n_sigma_vec = size(Chi, 2);
    % Propagated estimate
  x_bar = zeros(dim,1);
  for kk = 1:n_sigma_vec
    x_bar = x_bar + w_m(kk)*Chi(:,kk);
  end
    % Estimate covariance
  AT = [Chi(:,2:n_sigma_vec) Sr_Rv];
  for kk = 2:n_sigma_vec
    AT(:,kk-1) = sr_w_c(kk)*(Chi(:,kk) - x_bar);
  end
  [~, S_bar] = mth_qr(AT');
  S_bar = mth_chol_upd(S_bar, sr_w_c(1), Chi(:,1) - x_bar);
  S_bar = S_bar';
