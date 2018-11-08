function Chi = est_ut_sigma_vec(phat, SigmaPhat, alpha, kappa)
% EST_UT_SIGMA_VEC computes the unscented transform of an estimate
% given its covariance and "tuning" parameters.
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
%
% Return:
%   [nx(2n+1)] matrix of sigma vectors
%
% Kurt Motekew   2018/11/07
%

  n = size(phat,1);
  lambda = alpha*alpha*(n + kappa) - n;
  SigmaPscaled = sqrtm((n + lambda)*SigmaPhat);

  Chi = zeros(n, 2*n + 1);
  Chi(:,1) = phat;
  for ii = 2:(n+1)
    Chi(:,ii) = phat + SigmaPscaled(:,ii-1);
    Chi(:,ii+n) = phat - SigmaPscaled(:,ii-1);
  end

