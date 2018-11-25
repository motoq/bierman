function [x_bar, P_bar] = est_pred_ukf(Chi, w_m, w_c, Rv)
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
%   Chi   Sigma vectors to be condensed into a new estimate, [mxn]
%         where m is the dimension of the state vector and n is the
%         number of sigma vectors.
%   w_m   Estimate weighting factors, [1xn]
%   w_c   Covariance weighting factors, [1xn]
%
% Return:
%   x_bar  Updated estimate based on propagated sigma vectors, [mx1]
%   P_bar  Updated estimate covariance based on propagated sigma vectors, [mxm]
%
% Kurt Motekew   2018/11/14
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
  P_bar = zeros(dim);
  for kk = 1:n_sigma_vec
    chi_minus_xbar = Chi(:,kk) - x_bar;
    P_bar = P_bar + w_c(kk)*(chi_minus_xbar*chi_minus_xbar');
  end
    % Include process noise
  P_bar = P_bar + Rv;
