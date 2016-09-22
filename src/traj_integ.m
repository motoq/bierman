function [t, x] = traj_integ(t0, dt, tf, x0)
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
%   x0     Initial position and velocity, [6x1]
%
% Return
%   t   Array of times associated with each position, [1xN]
%   x   Position and velocity as a function of time, [6xN]
%
% Kurt Motekew   2014/11/06
%
  t = t0:dt:tf;
  ntimes = size(t,2);
  x = zeros(6,ntimes);
  x(1:6,1) = x0;
  for ii = 2:ntimes
    [~, x(:,ii)] = mth_rk4(@traj_dxdt, t(ii-1), dt, x(:,ii-1));
      % Below "ground" - trim and return
    if x(3,ii) < 0
      t = t(1:ii);
      x = x(:,1:ii);
      break;
    end
  end
  
