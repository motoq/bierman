function [newt, newx] = mth_rk4(f, tt, dt, xx)
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
%   f    Function to be integrated, f(tt, xx)
%   tt   Current time
%   dt   Time to integrate by
%   xx   Current state
%
% Return
%  newt   New time, tt + dt
%  newx   Updated state vector
%

  xd = f(tt, xx);
  xa = dt*xd;
   x = xx + 0.5*xa;
   t = tt + .5*dt;
    %
  xd = f(t, x);
   q = dt*xd;
   x = xx + .5*q;
  xa = xa + q + q;
    %
  xd = f(t, x);
   q = dt*xd;
   x = xx + q;
  xa = xa + q + q;
  newt = tt + dt;
    %
  xd = f(newt, x);
  newx = xx + (xa + dt*xd)/6.0;
