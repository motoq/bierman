function dx = traj_dxdt(t, x)
% TRAJ_DXDT computes the differential equation of a simple ballistic
% trajectory.  Gravity and the ballistic coefficient are read by calling
% external functions.
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
%   t      Time at which to evaluate derivative
%   x      Current state vector [6x1]
%            x[1] = pos_x
%            x[2] = pos_y
%            x[3] = pos_z
%            x[4] = vel_x
%            x[5] = vel_y
%            x[6] = vel_z
% 
% Return
%   dx     Derivatives to be integrated [6x1]
%            dx[1] = vel_x = x[4]
%            dx[2] = vel_y = x[5]
%            dx[3] = vel_z = x[6]
%            dx[4] = acc_x = b*vx + gx = b*x[4] + gx
%            dx[5] = acc_y = b*vy + gy = b*x[5] + gy
%            dx[6] = acc_z = b*vz + gz = b*x[6] + gz
%
% Kurt Motekew   2014/11/16
%

  b = traj_bcoeff();
  g = traj_gravt();

  v = x(4:6,1);

  dx(4:6,1) = b*v + g;
  dx(1:3,1) = v;

