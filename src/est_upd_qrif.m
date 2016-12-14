function [x, R, qty] = est_upd_qrif(R, qty, A, z, SqrtW)
% EST_UPD_QRIF Updates the apriori estimate and covariance via a QR IF
% (QR information filter) given a new set of observations.
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
%   R         R portion of QR decomposition performed on a priori A matrix,
%             [NxN]
%   qty       A priori Q*z, with Q portion of QR decomposition performed on
%             a priori A matrix, [Nx1]
%   A         A priori partial of obs w.r.t. current estimate, [MxN]
%   z         Observation, [Mx1]
%   SqrtW     Square root of observation weighting matrix, [MxM]
%
% Return:
%   x        Updated estimated, [Nx1]
%   R        Updated R, upper triangular, [NxN]
%   qty      Updated qty, [Nx1]
%
% Kurt Motekew   2016/08/11
%

  A = [R ; SqrtW*A];
  y  = [qty ; SqrtW*z];

  [Q, R] = mth_qr(A);
  qty = Q'*y;
  x = mth_trisol(R, qty);
