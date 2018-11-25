function [P] = mth_ud2p(UD)
% MTH_UD2P converts the combined U-D matrix UD to P = UDU'
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
%   UD  Combined unit diagonal upper triangular matrix U and diagonal matrix
%       D.  D is along the diagonal since U's diagonal would be 1's.  [NxN]
% Return:
%   P   UDU', [NxN]
%
% Author:  Kurt Motekew    20160803
% 
  n = size(UD,2);

  U = eye(n);
  D = zeros(n);
  for ii = 1:n
    D(ii,ii) = UD(ii,ii);
    for jj = (ii+1):n
      U(ii,jj) = UD(ii,jj);
    end
  end
  P = U*D*U';
