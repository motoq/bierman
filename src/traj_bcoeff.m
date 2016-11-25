function b = traj_bcoeff
% TRAJ_COEFF returns the ballistic "coefficient" used in this simulation
%
%-----------------------------------------------------------------------
% Copyright 2014 Kurt Motekew
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%-----------------------------------------------------------------------
%
% Return
%  b   Coefficient multiplied by velocity to result in acceleration,
%      should be negative for the acceleration to oppose the direction
%      of flight
%

  b = 0;%-.05;
