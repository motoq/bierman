function A = drng_dloc(pos, loc)
% DRNG_DLOC computes the partials of a range measurement w.r.t. a location
% in Cartesian coordinates given a the position of a range only tracker and
% a location
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
%   to the tracked object location [1x3]
%
% Kurt Motekew   2014/10/25
%
  srng = loc - pos;
  shat = srng/norm(srng);
  A = shat';
