function [Chi, w_m, w_c] = est_ut_sigma_vec(phat, SigmaPhat, alpha, kappa, beta)
% EST_UT_SIGMA_VEC computes the unscented transform of an estimate
% given its covariance and "tuning" parameters.  Also returns sigma vector
% and covariance weightings
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
%   phat       Estimate vector, [nx1]
%   SigmaPhat  1-sigma covariance of phat
%   alpha      Scaling parameter, small value typically 1e-4 <= alpha <=1
%   kappa      Secondary scaling parameter, often zero
%   beta       Distribution factor, 2 is optimal for Gaussian
%
% Return:
%   Chi  [nx(2n+1)] matrix of sigma vectors
%   w_m  Sigma vector weighting, [1x(2n+1)]
%   w_c  Sigma Covariance weighting, [1x(2n+1)]
%
% Kurt Motekew   2018/11/07
%
% Ref:  Rudolph van der Merwe & Eric A. Wan, "The Unscented Kalman Filter
%       for Nonlinear Estimation"
%

  n = size(phat,1);
  lambda = alpha*alpha*(n + kappa) - n;
  ScaledS = mth_sqrtm((n + lambda)*SigmaPhat);

    % Sigma vector
  Chi = zeros(n, 2*n + 1);
  Chi(:,1) = phat;
  for ii = 2:(n+1)
    Chi(:,ii) = phat + ScaledS(:,ii-1);
    Chi(:,ii+n) = phat - ScaledS(:,ii-1);
  end
  
    % Weighting
  w_0_m = lambda/(n + lambda);
  w_0_c = w_0_m - alpha*alpha + 1 + beta;
  w_i_m = 1/(2*(n + lambda));
  w_i_c = w_i_m;
    % Return as an array for future scalability and ease of use
  w_m(1,2:(2*n + 1)) = w_i_m;
  w_m(1,1)   = w_0_m;
  w_c(1,2:(2*n + 1)) = w_i_c;
  w_c(1,1)   = w_0_c;
