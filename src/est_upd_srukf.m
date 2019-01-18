function [x_hat, L_hat] = est_upd_srukf(x_bar, L_bar, Chi, w_m, sr_w_c,...
                                                           Y, y, Sr_Rn, nc)
% EST_UPD_SRUKF Given the current estimate and covariance, update with
% sigma vector based parameters and computed observations using a square root
% formulation of the unscented Kalman filter.  The input estimate may be an
% augmented state vector with consider parameters (this must be reflected in
% the input estimate covariance).  Consider parameters and the consider
% portion of the covariance (square root) will also not be altered.
%
%-----------------------------------------------------------------------
% Copyright 2019 Kurt Motekew
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%-----------------------------------------------------------------------
%
% Inputs:
%   x_bar    Current estimate [mX1]
%   L_bar    Estimate covariance, square root, LOWER triangular such
%            that P = L*L' [mXm]
%   Chi      Sigma vectors (used to form x_bar and L_bar), [mXn]
%            where n is the number of sigma vectors.
%   w_m      Estimate weighting factors, [1Xn]
%   sr_w_c   Square root of covariance weighting factors, [1Xn]
%   Y        Sigma vector based computed observations, [num_obs X n]
%   y        Observations, [num_obx X 1]
%   Sr_Rn    Square root of observation noise, [mxm]
%   nc       Number of parameters to consider but not update.  The
%            last nc entries in x_bar.  If zero, then you probably don't
%            know what you are doing and are trying to hack it with
%            process noise...
%
% Return:
%   x_hat  State estimate update based on observations
%   L_hat  Updated estimate covariance square root, lower triangular, such that
%          P = L_hat*L_hat', [mxm]
%
% Kurt Motekew   2019/01/15
%

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
    % Compute Kalman gain
  [~, S_Y_bar] = mth_qr(AT');
  R_Y_bar = mth_triinv(S_Y_bar);
  K = SigmaXY*(R_Y_bar*R_Y_bar');
    % Apply to estimate covariance square root via successive
    % rank one downdates
  U = K*S_Y_bar;
  n = size(U, 2);
  L_hat = L_bar;
  for kk = 1:n
    L_hat = mth_chol_upd_l(L_hat, -1, U(:,kk));
  end
    % If consider parameters are present, reset that portion of
    % the covariance square root via successive rank one updates...
  if (nargin == 9)  &&  (nc > 0)
    np = dim-nc;
    Kc = [zeros(np,n_obs) ; K((np+1):dim,:)];
    U2 = Kc*S_Y_bar;
    for kk = 1:n
      L_hat = mth_chol_upd_l(L_hat, 1, U2(:,kk));
    end
      % ...and modify the Kalman gain to not affect the consider parameters
    K((np+1):dim,:) = zeros(nc,n_obs);
    x_hat = x_bar + K*(y - y_bar);
  else
      % ...otherwise just knock it out,
    x_hat = x_bar + K*(y - y_bar);
  end
 
