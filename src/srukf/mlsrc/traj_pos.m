function pos = traj_pos(dt, pos0, vel0)
% TRAJ_POS computes position as a function of time given initial position and
% velocity assuming only a constant gravitational acceleration via
% traj_gravt influences
%
%-----------------------------------------------------------------------
% Copyright 2014 Kurt Motekew
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%-----------------------------------------------------------------------
%
%   x = x0 + v0*t + .5*a*t^2
%
% Inputs:
%   dt    Time from epoch for which to compute a new position
%   pos0  Initial position, [3x1]
%   vel0  Initial velocity, [3x1]
%
% Return:
%   Updated position, [3x1]
%
  acc = traj_gravt();
  pos = pos0 + vel0*dt + 0.5*acc*dt*dt;
