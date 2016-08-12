function [A] = mth_householder_tri(A)
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
%   A  Triangularized, [MxN] matrix.
%
% Author:  Kurt Motekew    20160809
% 

  [m, n] = size(A);

  for ii = 1:n
      % Already processed last row - can't go farther
    if ii > m
      break;
    end
    A(ii:m,ii:n) = mth_hh_tri_col(A(ii:m,ii:n));
  end

