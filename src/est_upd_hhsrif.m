function [xhat, Rhat, zhat, e] = est_upd_hhsrif(R0, z0, A, z, SqrtW)
% EST_UPD_HHSRIF Updates the apriori estimate and covariance via a SRIF
% using Householder triangularization given a new set of observations
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
%   R0        A priori partials component of triangularized information array,
%             [NxN]
%   z0        A priori obs portion of triangularized info arry, [Nx1]
%   A         Partial of obs w.r.t. current estimate, [MxN]
%   z         Observation vector, [Mx1]
%   SqrtW     Square root of observation weighting matrix, [MxM]
%
% Return:
%   xhat     Updated estimated, [Nx1]
%   Rhat     Updated R0, upper triangular, [NxN]
%   zhat     Updated z0, [Nx1]
%   e        Residual, [(M-N)x1]
%
% Kurt Motekew   2016/08/11
%
%
% Ref:  G. J. Bierman, Factorization Methods for
%       Discrete Sequential Estimation, Dover Publications, Inc.,
%       Mineola, NY, 1977, pp. 71-72.
%

    % Normalize when forming info array, Covariance obs = I
  R_z = [ R0 z0 ; SqrtW*A SqrtW*z ];
  R_z_hat = mth_householder_tri(R_z, 1);

  n = size(z0,1);
  m = size(z,1);
  Rhat = R_z_hat(1:n,1:n);
  zhat = R_z_hat(1:n,(n+1));
  e = R_z_hat((n+1),(n+1):m);
  
  xhat = mth_trisol(Rhat, zhat);
