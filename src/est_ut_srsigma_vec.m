function Chi = est_ut_srsigma_vec(phat, Shat, alpha, kappa)
% EST_UT_SRSIGMA_VEC computes the unscented transform of an estimate
% given the square root of its covariance and "tuning" parameters.
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
%   phat   Estimate vector, [nx1]
%   Shat   Square root of 1-sigma covariance of phat, such that P = S*S' [nxn]
%   alpha  Scaling parameter, small value typically 1e-4 <= alpha <=1
%   kappa  Secondary scaling parameter, often zero
%
% Return:
%   Chi  [nx(2n+1)] matrix of sigma vectors
%
% Kurt Motekew   2018/12/25
%
% Ref:  Rudolph van der Merwe & Eric A. Wan, "The Unscented Kalman Filter
%       for Nonlinear Estimation"
%

  n = size(phat,1);
  lambda = alpha*alpha*(n + kappa) - n;
  ScaledS = sqrt(n + lambda)*Shat;

    % Sigma vector
  Chi = zeros(n, 2*n + 1);
  Chi(:,1) = phat;
  for ii = 2:(n+1)
    Chi(:,ii) = phat + ScaledS(:,ii-1);
    Chi(:,ii+n) = phat - ScaledS(:,ii-1);
  end
