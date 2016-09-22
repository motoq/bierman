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
x0 = [0.35 0.25 .25 0.2  0.2  1]';

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

  % Create an ideal trajectory.  This is the model used by the filter
  % without the influence of observations.
[~, x_ideal] = traj_ideal(t0, dt, tf, x0);
[t, x_integ] = traj_integ(t0, dt, tf, x0);
  % number of observed state vectors over time of interest or until impact
nobs = size(t,2);

  % Plot geometry
traj_plot(x_ideal, x_integ, tkrs, blen);
title('Boxed Tracker Geometry');
view([70 20]);

%
% Create simulated range measurements and plot
%

z = zeros(nmeas,nobs);
for ii = 1:nobs
  for jj = 1:nmeas
    svec = x_integ(1:3,ii) - tkrs(:,jj);
    z(jj,ii) = norm(svec) + srng*randn;
  end
end
figure;
mesh(z);
xlabel('Observation');
ylabel('Tracker');
title('Measured Range to Tracked Object');



