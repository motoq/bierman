function vel = traj_vel(dt, vel0)
% TRAJ_VEL computes velocity as a function of time given initial velocity
% assuming only gravitational acceleration via traj_gravt influences
%
%-----------------------------------------------------------------------
% Copyright 2014 Kurt Motekew
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%-----------------------------------------------------------------------
%
%   v = v0 +  a*t
%
% Inputs:
%   dt     Time from epoch for which to compute the new velocity
%   vel0   Initial velocity
%
% Return:
%   Updated velocity
%
  acc = traj_gravt();
  vel = vel0 + acc*dt;
