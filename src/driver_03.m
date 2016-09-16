%
% Copyright 2016 Kurt Motekew
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%

% Beginning of dynamic scenario - smacking a wiffle ball in a room trajectory

close all;
clear;

  %
  % SETUP
  %

  % Tracked object initial state and trajectory
rho0 = [0.35 0.25 .25]';
vel0 =  [0.2  0.2  1]';

  % Tracker locations
blen = 1;
tkrs = [
         0       0       blen ;
         blen    0       blen ;
         blen    blen    blen ;
         0       blen    blen
      ]';
  % Number of trackers = number of measurements per observation set
  % Observation sets are obs at a single instant in time
nmeas = size(tkrs,2);
  % Tracker accuracy
srng = .001;

t0 = 0;                                          % Scenario start time
tf = 3;                                          % Scenario stop time
dt = .01;                                        % Data rate

  % Create an ideal trajectory, as assumed by the filter,
  % and the actual trajectory for which measurements will be based
[~, rhos_ideal, vels_ideal] = traj_ideal(t0, dt, tf, rho0, vel0);
[t, rhos, vels] = traj_integ(t0, dt, tf, rho0, vel0);
  % number of observed state vectors over time of interest or until impact
nobs = size(t,2);

  % Plot geometry
traj_plot(rhos_ideal, rhos, tkrs, blen);
title('Boxed Tracker Geometry');
view([70 20]);


