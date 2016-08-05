function ApTWAp = box_infom(r, rho, W)
% BOX_INFOM computes the information matrix associated with a postulated
% tracked object location and range only trackers.
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
%   r     A [3xN] matrix of N tracker positions
%   rho   The postulated position of the tracked object, [3x1]
%   W     Range uncertianty weighting matrix, 1/sigma_range^2, scalar
%
% Return:
%   ApTWAp   [3x3] information matrix built from the partials of
%            range measurements w.r.r. the tracked object's unknown
%            coordinates
%
% Kurt Motekew   2014/10/25
%

  nmeas = size(r,2);
  ApTWAp = zeros(3);
  for ii = 1:nmeas
    Api = drng_dloc(r(:,ii), rho);
    ApTWAp = ApTWAp + Api'*W*Api;
  end

