function [t, x] = mth_rk4(f, t0, dt, x0)
% MTH_RK4 performs Runge-Kutta Integration
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
%   f    Function to be integrated, f(t0, x0)
%   t0   Current time
%   dt   Time to integrate by
%   x0   Current state
%
% Return
%  t   New time, t0 + dt
%  x   Updated state vector
%

  xd = f(t0, x0);
  xa = dt*xd;
   x = x0 + 0.5*xa;
   t = t0 + .5*dt;
    %
  xd = f(t, x);
   q = dt*xd;
   x = x0 + .5*q;
  xa = xa + q + q;
    %
  xd = f(t, x);
   q = dt*xd;
   x = x0 + q;
  xa = xa + q + q;
   t = t0 + dt;
    %
  xd = f(t, x);
  x = x0 + (xa + dt*xd)/6.0;
