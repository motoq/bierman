function A = est_drpnt_dloc(pos, loc)
% EST_DRPNT_DLOC computes the partials of a range and pointing measurements
% w.r.t. a location in Cartesian coordinates given a the position of a
% the tracker and the location being tracked.
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
%   Jacobian relating the partials of the slant range distance and pointing
%   measurements from the tracker to the tracked object location [3x3]
%   The pointing measurements are the x and y components of a unit pointing
%   vector from the tracker.
%
% Kurt Motekew   2016/12/11
%
    % Pointing and range
  svec = loc - pos;
  smag = norm(svec);
  sinv = 1/smag;
  shat = sinv*svec;

  A = zeros(3);
  A(1,:) = shat';                                % First row is range
  Prj_scaled = sinv*(eye(3) - shat*shat');       % Scaled projection
  A(2:3,:) = Prj_scaled(1:2,:);                  % 2nd and 3rd are pointing
