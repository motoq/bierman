function [t, rhos, vels] = traj_integ(t0, dt, tf, rho0, vel0)
% TRAJ_INTEG computes a ballistic trajectory with drag effects given initial
% position and velocity using the differential equations supplied by
% traj_dxdt using RK4 numeric integration:
%
%     ddx/ddt = b*dx/dt + g
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
  rhos(:,1) = rho0;
  vels(:,1) = vel0;
  rho = rho0;
  vel = vel0;
  for ii = 2:ntimes
    [~, xx] = mth_rk4(@traj_dxdt, t(ii), dt, [ rho' vel' ]');
    rho = xx(1:3, 1);
    vel = xx(4:6, 1);
    if rho(3,1) < 0
      ntimes = ii-1;
      rhos = rhos(:,1:ntimes);
      vels = vels(:,1:ntimes);
      t = t(1:ntimes);
      break;
    else
      rhos(:,ii) = rho;
      vels(:,ii) = vel;
    end
  end
  
