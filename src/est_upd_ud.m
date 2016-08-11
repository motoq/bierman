function [x, U, D] = est_upd_ud(x, U, D, Ap, z, r)
% EST_UPD_UD Updates the apriori estimate and covariance via U-D sequential
% estimation given a single new observation.  The covariance P is input in
% the decomposed form such that P = UDU' where U is a unit upper triangular
% matrix and D is a diagonal matrix.
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
%   x         A priori estimate, [Nx1]
%   U         A priori U matrix, [NxN]
%   D         A priori D matrix, [NxN]
%   Ap        Partial of obs w.r.t. x, [1xN]
%   z         Observation residual, scalar
%   r         Observation variance, sigma^2, scalar
%
% Return:
%   x   Updated estimated, [Nx1]
%   U   Updated U, [NxN]
%   D   Updated D, [NxN]
%
% Kurt Motekew   2016/08/03
%
    %
    % Ref:  G. J. Bierman, Factorization Methods for
    %       Discrete Sequential Estimation, Dover Publications, Inc.,
    %       Mineola, NY, 1977.
    %
  n = size(x,1);                                 % Number of solve for

  a = U'*Ap';
  b = D*U'*a;
  
  alpha = r + b(1)*a(1);
  gamma = 1/alpha;
  D(1,1) = r*gamma*D(1,1);

  for jj = 2:n
    beta = alpha;
    alpha = alpha + b(jj)*a(jj);
    lambda = -a(jj)*gamma;
    gamma = 1/alpha;
    D(jj,jj) = beta*gamma*D(jj,jj);
    for ii = 1:(jj-1)
      beta = U(ii,jj);
      U(ii,jj) = beta + b(ii)*lambda;
      b(ii) = b(ii) + b(jj)*beta;
    end
  end
  x = x + b*z*gamma;
