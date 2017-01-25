function A = est_drng_dpos(pos, loc)
% EST_DRNG_DPOS computes the partials of a range measurement w.r.t. the
% range tracker position in Cartesian coordinates.
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
%   pos   Tracker location, [3x1]
%   loc   Location of tracked object, [3x1]
%
% Return:
%   Jacobian relating the partials of the slant range distance
%   to the tracker location [1x3]
%
% Kurt Motekew   2017/01/02
%
  srng = pos - loc;
  shat = srng/norm(srng);
  A = shat';
