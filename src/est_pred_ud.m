function [xbar, Ubar, Dbar] = est_pred_ud(xhat, Uhat, Dhat, Phi, Q, G)
% EST_PRE_UD Updates a prior estimate's U-D covariance given a state
% transition matrix and process noise matrix.  The decomposed form is
% such that
%   P = UDU'
% where U is a unit upper triangular matrix and D is a diagonal matrix.
% Note that the process noise matrix is diagonal.  Consult Bierman,
% p. 125, w.r.t. whitening a non-diagonal Q matrix.
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
%   xhat      A priori estimate, [nx1]
%   Uhat      A priori U matrix, [nxn]
%   Dhat      A priori D matrix, [nxn]
%   Phi       State transition matrix, [nxn]
%   Q         Process noise matrix, DIAGONAL of square mxm matrix, [1xm] 
%   G         Partials of estimate w.r.t. process noise, dx/dq, [nxm]
%
% Return:
%   xbar   Updated estimated, [nx1]
%   Ubar   Updated U, [nxn]
%   Dbar   Updated D, [nxn]
%
% Kurt Motekew   2016/12/14
%
%
% Ref:  G. J. Bierman, Factorization Methods for
%       Discrete Sequential Estimation, Dover Publications, Inc.,
%       Mineola, NY, 1977
%       
%       Byron D. Tapley, Bob E. Schutz, George H. Born, Statistical
%       Orbit Determination, Elsevier Academic Press, 2004.
%         The algorithm from page 345 was invaluable in deciphering
%         a form that could be represented in terms of matrix algebra.
%

  xbar = Phi*xhat;

    % Get matrix sizes for looping
  n = size(Uhat,1);                         % Number of solve for
  m = size(Q,2);                            % Number of consider parameters
  N = n + m;

  A = [ Phi*Uhat G]';                       % [Nxn], a_j = A(:,j)
  Dtmp((n+1):N,(n+1):N) = diag(Q);          % [NxN]
  Dtmp(1:n,1:n) = Dhat;

    % Could reuse Uhat and Dhat
  Ubar = eye(n);
  Dbar = zeros(n);
  for ii = n:-1:1
    c = Dtmp*A(:,ii);
    Dbar(ii,ii) = A(:,ii)'*c;
    d = c/Dbar(ii,ii);
    for jj = 1:n
      Ubar(jj,ii) = A(:,jj)'*d;
      A(:,jj) = A(:,jj) - Ubar(jj,ii)*A(:,ii);
    end
  end

    % Cheating method - but good for verification
    % This is what we are doing in a more recursive manner
%  P = U*D*U';
%  P = Phi*P*Phi' + G*diag(Q)*G';
%  [U, D] = mth_udut2(P);
