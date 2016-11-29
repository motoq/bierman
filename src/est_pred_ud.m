function [x, U, D] = est_pred_ud(x, U, D, Phi, Q, G)
% EST_PRE_UD Updates a prior estimate's U-D covariance given a state transition
% matrix and process noise matrix.  The decomposed form is such that
% P = UDU' where U is a unit upper triangular matrix and D is a diagonal matrix.
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
%   x         A priori estimate, [nx1]
%   U         A priori U matrix, [nxn]
%   D         A priori D matrix, [nxn]
%   Phi       State transition matrix, [nxn]
%   Q         Process noise matrix, Diagonal of square Nw matrix, [1xNw] 
%   G         Partials of estimate w.r.t. process noise, dx/dq, [nxNw]
%
% Return:
%   x   Updated estimated, [Nx1]
%   U   Updated U, [NxN]
%   D   Updated D, [NxN]
%
% Kurt Motekew   2016/11/22
%
%
% Ref:  G. J. Bierman, Factorization Methods for
%       Discrete Sequential Estimation, Dover Publications, Inc.,
%       Mineola, NY, 1977, pp. 132-133
%

    % Linear state propagation
  x = Phi*x;

    % Cheat for now - need to implement recursive form
  P = U*D*U';
  P = Phi*P*Phi' + G*diag(Q)*G';
  [U, D] = mth_udut2(P);

    % Get matrix sizes for looping
%  n = size(Phi,1);
%  Nw = size(Q,2); 
%  N = n + Nw;

%  v = zeros(1,n);
  
