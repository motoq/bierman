function [R, b] = est_pred_hhsrif(R, b, PhiInv, Rw, G)
% EST_PRE_HHSRIf Updates a prior estimate's square root information matrix
% in time give a state transition matrix.  Also incorporates proces noise.
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
%   R         a priori estimate information matrix square root, [nxn]
%   b         a priori observation, [nx1]
%   PhiInv    Inverse of state transition matrix, [nxn]
%   Rw        Inverse of square root of process noise covariance, [mxm]
%   G         Partials of estimate w.r.t. process noise, dx/dq, [nxm]
%
% Return:
%   R      Updated estimate information matrix square root, [nxn]
%   b      Updated observation, [nx1]
%
% Kurt Motekew   2016/12/18
%
%
% Ref:  G. J. Bierman, Factorization Methods for
%       Discrete Sequential Estimation, Dover Publications, Inc.,
%       Mineola, NY, 1977, pp. 115-122.
%

    % Get matrix sizes for indexing
  n = size(R,2);                            % Number of solve for
  m = size(Rw,1);                           % Number of consider parameters

  R = R*PhiInv;
  A = [ Rw           zeros(m,n)  zeros(m,1) ;
       -R*G          R           b           ];
    % Perform m transformations
  Ahat = mth_householder_tri(A, n+1);

    % Locate only R and b for now
  r1 = m+1;
  c1 = m+1;
  R = Ahat(r1:(r1+n-1),c1:(c1+n-1));
  b = Ahat(r1:(r1+n-1),c1+n);

