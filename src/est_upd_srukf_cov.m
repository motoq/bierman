function L_hat = est_upd_srukf_cov(x_bar, L_bar, Chi, w_m, sr_w_c, Y, Sr_Rn, nc)
% EST_UPD_SRUKF_COV Given a postulated estimate and covariance, update with
% sigma vector based parameters and computed observations using a square root
% formulation of the unscented Kalman filter.  The input estimate may be an
% augmented state vector with consider parameters (this must be reflected in
% the input estimate covariance).  The consider parameter portion of the
% covariance (square root) will not be altered.
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
%   x_bar    Postulated estimate [mX1]
%   L_bar    Estimate covariance, square root, LOWER triangular such
%            that P = L*L' [mXm]
%   Chi      Sigma vectors (used to form x_bar and L_bar), [mXn]
%            where n is the number of sigma vectors.
%   w_m      Estimate weighting factors, [1Xn]
%   sr_w_c   Square root of covariance weighting factors, [1Xn]
%   Y        Sigma vector based computed observations, [num_obs X n]
%   Sr_Rn    Square root of observation noise, [mxm]
%   nc       Number of parameters to consider but not update.  The
%            last nc entries in x_bar.
%
% Return:
%   L_hat  Updated estimate covariance square root, lower triangular, such
%          that P = L_hat*L_hat', [mxm]
%
% Kurt Motekew   2019/01/27
%
% Ref:  Jeroen L. Geeraert & Jay W. McMahon, "Square-Root Unscented
%       Schmidt-Kalman Filter"
%
% Note this implementation of the square root UKF assumes the covariance
% weighting factors are all positive (see the input sr_w_c).
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
    % Apply to estimate covariance square root via successive
    % rank one downdates
  U = l_mth_kst(SigmaXY, S_Y_bar);
  n = size(U, 2);
  L_hat = L_bar;
  for kk = 1:n
    L_hat = mth_chol_upd_l(L_hat, -1, U(:,kk));
  end
    % If consider parameters are present, reset that portion of
    % the covariance square root via successive rank one updates...
  if (nargin == 8)  &&  (nc > 0)
    np = dim-nc;
    U(1:np,1:n_obs) = 0;
    for kk = 1:n
      L_hat = mth_chol_upd_l(L_hat, 1, U(:,kk));
    end
  end
end
 

function Kst = l_mth_kst(P, S)
% L_MTH_KST Uses forward substitution to solve for K*S' given P and S such
% that K*S'*S = P and S is upper triangular.
%
  [m, n] = size(P);
    % Solve for K*S'
  Kst = zeros(m,n);
  for ii = 1:m
    Kst(ii,1) = P(ii,1)/S(1,1);
    for jj = 2:n
      s = P(ii,jj);
      for kk = jj:-1:2
        s = s - Kst(ii,kk-1)*S(kk-1,jj);
      end
      Kst(ii,jj) = s/S(jj,jj);
    end
  end
end
