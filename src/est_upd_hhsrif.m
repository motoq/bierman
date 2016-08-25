function [phat, R, z] = est_upd_hhsrif(phat0, R0, z0,  Ap, r, sw)
% EST_UPD_HHSRIF Updates the apriori estimate and covariance via a SRIF
% using Householder triangularization given a single new observation
%
%-----------------------------------------------------------------------
% Copyright 2016 Kurt Motekew
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%-----------------------------------------------------------------------
%
% Inputs:
%   phat0     A priori estimate, [Nx1]
%   R0        A priori partials component of triangularized information array,
%             [NxN]
%   z         A priori obs portion of triangularized info arry, [Nx1]
%   Ap        Partial of obs w.r.t. phat0, [1xN]
%   r         Observation residual, scalar
%   sw        Inverse of observation uncertainty, 1/sigma, scalar
%
% Return:
%   phat     Updated estimated, [Nx1]
%   R        Updated R0, [NxN]
%   z        Updated z0, [Nx1]
%
% Kurt Motekew   2016/08/11
%
%
% Ref:  G. J. Bierman, Factorization Methods for
%       Discrete Sequential Estimation, Dover Publications, Inc.,
%       Mineola, NY, 1977, pp. 71-72.
%

    % Normalize when forming info array, Covariance obs = I
  R_z = [ R0 z0 ; sw*Ap sw*r ];
  R_z_hat = mth_householder_tri(R_z, 1);

  n = size(phat0,1);
  R = R_z_hat(1:n,1:n);
  z = R_z_hat(1:n,(n+1));
  
  dp = mth_trisol(R, z);
  phat = phat0 + dp;
