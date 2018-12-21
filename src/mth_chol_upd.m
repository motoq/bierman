function S = mth_chol_upd(S, c, a)
% MTH_CHOL_UPD Given S = sqrt(P) such that P = S*S', compute S2 such
% that P2 = P + c*a*a' where 'c' is a scalar, 'a' is a vector, and
% P2 = S2*S2'.
%
%-----------------------------------------------------------------------
% Copyright 2018 Kurt Motekew
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%-----------------------------------------------------------------------
%
% Inputs:
%   S   Upper triangular square root of an [NxN] matrix P such that
%       P = S*S'
%   c   Scalar
%   a   Vector [Nx1]
% Return:
%   S   Updated upper triangular square root such that
%       S = sqrt(S*S' +  c*a*a')
%
% Author:  Kurt Motekew    20181220
% 

  [U, D] = mth_s2udut2(S);

  [U, D] = mth_udut2_upd(U, D, c, a);
  S = U*sqrt(D);
