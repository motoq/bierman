function [phat, R, qty] = est_upd_qrsrif(phat0, R0, qty0,  Ap, r, sw)
% EST_UPD_QRSRIF Updates the apriori estimate and covariance via a SRIF
% using QR decomposition given a single new observation
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

  A = [R0 ; sw*Ap];
  y  = [qty0 ; sw*r];

  %[Q, R] = qr(A, '0');
  [Q, R] = mth_qr(A);
  qty = Q'*y;
  dp = mth_trisol(R, qty);
  phat = phat0 + dp;
