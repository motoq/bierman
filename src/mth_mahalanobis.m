function md = mth_mahalanobis(p, phat, SigmaP)
% MTH_MAHALANOBIS computes the Mahalanobis distance given vectors of true
% and estimated values, and the associated estimate covariance.
%
%-----------------------------------------------------------------------
% Copyright 2016 Kurt Motekew
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%-----------------------------------------------------------------------
%
% Inputs
%   p        Vector of true values, [Nx1]
%   phat     Vector of estimated values, [Nx1]
%   SigmaP   Estimate covariance, [NxN]
%
% Kurt Motekew   2014
%

  r = p - phat;
  md = sqrt(r' * SigmaP^-1 *r);
