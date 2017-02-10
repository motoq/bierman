function [xhat, Rxhat, Rxyhat, zxhat, Ryhat, zyhat] =...
                       est_upd_hhsrif_bias(Rx, Rxy, zx, Ax, z, SqrtW,...
                                           Ry, zy, Ay)
% EST_UPD_HHSRIF_BIAS Updates the apriori estimate and covariance via a
% SRIF using Householder triangularization given a new set of observations
% and bias statistics
%
%-----------------------------------------------------------------------
% Copyright 2017 Kurt Motekew
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%-----------------------------------------------------------------------
%
% Inputs:
%   Rx        A priori triangularized information array, obs w.r.t. solve
%             for x, [NxN], N solve for
%   Rxy       A priori cross correlation partials, resulting from previous
%             calls to this function, [NxNy]
%   zx        A priori compressed obs, [Nx1]
%   Ax        Partial of obs w.r.t. current estimate, [MxN]
%   z         Observation vector, [Mx1]
%   SqrtW     Square root of observation weighting matrix, [MxM]
%   Ry        A priori partials of triangularized information array, obs
%             w.r.t. bias y, [NyxNy]
%   zy        A priori compressed bias, [Nx1]
%   Ay        Partial of obs w.r.t. bias, [MxNy]
%
% Return:
%   xhat     Updated estimated, [Nx1]
%   Rxhat    Updated Rx, upper triangular, [NxN]
%   Rxyhat   Updated Rxy, [NxNy]
%   zxhat    Updated z, [Nx1]
%   Ryhat    Updated Ry, NyxNy]
%   zyhat    pdated y
%   e        Sum of square of residuals, scalar
%
% Kurt Motekew   2017/01/22
%
%
% Ref:  G. J. Bierman, Factorization Methods for
%       Discrete Sequential Estimation, Dover Publications, Inc.,
%       Mineola, NY, 1977, pp. 71-72.
%

    % Normalize when forming info array, Covariance obs = I
  [nx, ny] = size(Rxy);

  Rxyz = [ Rx Rxy zx ;  SqrtW*Ax SqrtW*Ay SqrtW*z ];
  Rxyz_hat = mth_householder_tri(Rxyz, ny+1);
  rows = 1:nx;
  cols = 1:nx;
  Rxhat = Rxyz_hat(rows,cols);
  cols = (nx+1):(nx+ny);
  Rxyhat = Rxyz_hat(rows,cols);
  cols = nx + ny + 1;
  zxhat = Rxyz_hat(rows,cols);
  xhat   = mth_trisol(Rxhat, zxhat);

  m = size(z,1);
  rows = (1:m) + nx;
  cols = (1:ny) + nx;
  Ayhat = Rxyz_hat(rows,cols);
  cols = nx + ny + 1;
  zhat = Rxyz_hat(rows,cols);
  Ryz = [ Ry zy ; Ayhat zhat ];
  Ryz_hat = mth_householder_tri(Ryz, 1);
  Ryhat = Ryz_hat(1:ny,1:ny);
  zyhat = Ryz_hat(1:ny,ny+1);

  Rxhatinv = mth_triinv(Rxhat);
  Shat = -Rxhatinv*Rxyhat;
  yhat = mth_trisol(Ryhat, zyhat);
  xhat = xhat + Shat*yhat;
