function [t, x] = traj_ideal(t0, dt, tf, x0)
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
%   x0     Initial position and velocity, [6x1]
%    
% Return
%   t      Array of times associated with each position, [1xN]
%   x      Position and velocity vectors as a function of time, [3xN]
%
% Kurt Motekew   2014/11/06
%

  t = t0:dt:tf;
  ntimes = size(t,2);
  x = zeros(6,ntimes);
  for ii = 1:ntimes
    x(1:3,ii) = traj_pos(t(ii), x0(1:3), x0(4:6));         % Position update
    x(4:6,ii) = traj_vel(t(ii), x0(4:6));                  % Velocity update
      % Below the "ground" - done
    if x(3,ii) < 0
      ntimes = ii-1;
      x = x(:,1:ntimes);
      t = t(1:ntimes);
      break;
    end
  end
  
