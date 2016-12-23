function [TA] = mth_hh_tri_col(A)
% MTH_HH_TRI_COL Triangularizes the first column of an [MxN] matrix via a
% Householder Transformation.  Call this function once per column (working
% down row by row) to fully triangularize a matrix.
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
%   A   [MxN] Matrix
% Return:
%   TA  A with the first column Triangularized, [MxN] matrix.
%
% Author:  Kurt Motekew    20160809
%          Kurt Motekew    20161222  Vectorized some stuff but left and
%                                    commented out original for ease of
%                                    translation to other languages.
% 

  [m, n] = size(A);
  TA = zeros(m,n);

  u = zeros(1,m);

    % Diagonal
  %rssa = 0;
  %for ii = 1:m
  %  rssa = rssa + A(ii,1)*A(ii,1);
  %end
  rssa = sum(A(:,1).*A(:,1));
  rssa = sqrt(rssa);
  s = -sign(A(1,1))*rssa;
  TA(1,1) = s;

    % Donesky - column matrix - just compute the scalar along the "diagonal"
  if n == 1
    return;
  end

    % Reflect about u
  u(1) = A(1,1) - s;
  u(2:m) = A(2:m,1)';
  %for ii = 2:m
  %  u(ii) = A(ii,1);
  %end
  beta = 1/(s*u(1));

  for jj = 2:n
    a = zeros(1,m);
    ua = 0;
    for ii = 1:m
      ua = ua + u(ii)*A(ii,jj);
    end
    gamma = beta * ua;
    for ii = 1:m
      a(ii) = A(ii,jj) + gamma*u(ii);
    end
    TA(:,jj) = a;
  end

