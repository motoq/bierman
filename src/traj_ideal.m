function [t, rhos, vels] = traj_ideal(t0, dt, tf, rho0, vel0)
% TRAJ_IDEAL computes an ideal ballistic trajectory given initial position
% and velocity assuming no drag.  State vectors will be produced from the
% start until the end time, or until the z component of a state vector falls
% below zero.  See traj_pos and traj_vel for modeling assumptions.
%
%-----------------------------------------------------------------------
% Copyright 2014 Kurt Motekew
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%-----------------------------------------------------------------------
%
% Inputs
%   t0     Trajectory start time
%   dt     Time step at which to output state vectors
%   tf     Trajectory end time
%   rho0   Initial position, [3x1]
%   vel0   Initial velocity, [3x1]
%    
% Return
%   t      Array of times associated with each position, [1xN]
%   rhos   Positions as a function of time, [3xN]
%   vels   Velocities as a function of time, [3xN]
%
% Kurt Motekew   2014/11/06
%

  t = t0:dt:tf;
  ntimes = size(t,2);
  rhos = zeros(3,ntimes);
  vels = zeros(3,ntimes);
  for ii = 1:ntimes
    rhos(:,ii) = traj_pos(t(ii), rho0, vel0);
    vels(:,ii) = traj_vel(t(ii), vel0);
    if rhos(3,ii) < 0
      ntimes = ii-1;
      rhos = rhos(:,1:ntimes);
      vels = vels(:,1:ntimes);
      t = t(1:ntimes);
      break;
    end
  end
  
