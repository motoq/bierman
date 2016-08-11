function [TA] = mth_householder_tri(A)
% MTH_HOUSEHOLDER_TRI Triangularizes an [MxN] matrix via the Householder
% Transformation.
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
%   TA  Triangularized, [MxN] matrix.
%
% Author:  Kurt Motekew    20160809
% 

  [m, n] = size(A);
  TA = zeros(m,n);

  for ii = 1:n
    if ii > m
      break;
    end
    TA(ii:m,ii:n) = mth_hh_tri_col(A(ii:m,ii:n));
    A = TA;
  end

