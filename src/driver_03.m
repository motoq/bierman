%
% Copyright 2016 Kurt Motekew
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%

%
% Simple wiffle ball in a room trajectory
%
% This driver script is the first of this series of examples illustrating
% filtering of a dynamic scenario.  First, the "true" trajectory based on
% integrating a simple flat earth gravity model with a fixed drag coefficient
% is generated along with an ideal analytic model used by the filter that does
% not include drag is computed and plotted.  The range tracker locations are
% also included for reference.  The range uncertainty is then used to create
% simulated observations.  Various methods are then used to estimate the
% trajectory.
%
% 1) Batch Method (observation model only).  A WLS method is used to
%    process range values at each time step.  The current position estimate
%    and two prior estimates, along with the corresponding covariances, are
%    then used to estimate the velocity and velocity covariance associated
%    with the latest position estimate.
% 
% 2) 
%
% Kurt Motekew  2016/09/28
%

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
ntkrs = size(tkrs,2);
  % Tracker accuracy
srng = .001;
vrng = srng*srng;
W = vrng^-1;

t0 = 0;                                          % Scenario start time
tf = 3;                                          % Scenario stop time
dt = .01;                                        % Data rate

  % Create an ideal trajectory.  This is the model used by the filter
  % without the influence of observations.
[~, x_ideal] = traj_ideal(t0, dt, tf, x0);
  % 'True' trajectory
[t, x_true] = traj_integ(t0, dt, tf, x0);
  % number of measurement sets - total obs = nsets*ntkrs
nsets = size(t,2);

  % Plot geometry
traj_plot(x_ideal, x_true, tkrs, blen);
title('Boxed Tracker Geometry');
view([70 20]);

%
% Create simulated range measurements and plot
%

  % Vectors of measurement sets.  Columns at different times
z = zeros(ntkrs,nsets);
for ii = 1:nsets
  for jj = 1:ntkrs
    svec = x_true(1:3,ii) - tkrs(:,jj);
    z(jj,ii) = norm(svec) + srng*randn;
  end
end
figure;
mesh(z);
xlabel('Observation');
ylabel('Tracker');
title('Measured Range to Tracked Object');

  % Prime filters with batch derived state.  Use first three measurement
  % sets to generate three positions.  Then estimate velocity using the
  % three position estimates.  This creates an initial state for the third
  % time at which data is recorded.
ninit = 3;                                  % Number of initial batch estimates
x_0 = zeros(3, ninit);                      % Position vectors
P_0(3,3,ninit) = 0;                         % Position covariances
for ii = 1:ninit
  [x_0(:,ii), P_0(:,:,ii), ~] = box_locate(tkrs, z(:,ii)', W);
  %mth_mahalanobis(x_true(1:3,ii), x_0(1:3,ii), P_0(:,:,ii))
end
[vel, SigmaV] = traj_pos2vel(dt, x_0(:,1),  P_0(:,:,1),...
                                 x_0(:,2),  P_0(:,:,2),...
                                 x_0(:,3),  P_0(:,:,3));
%mth_mahalanobis(x_true(4:6,3), vel, SigmaV)

  % Number of filtered outputs will be the initial batch derived
  % estimate plus the remaining number of observations
nfilt = nsets - ninit + 1;
x = zeros(6,nfilt);
P = zeros(6,6,nfilt);
x(:,1) = [x_0(:,ninit) ; vel]; 
P(1:3,1:3,1) = P_0(:,:,3);
P(4:6,4:6,1) = SigmaV;

  % Batch method - observation model only
  % Don't use the first ninit observations.
for ii = 2:nfilt
  [x(1:3,ii), P(1:3,1:3,ii)] = box_locate(tkrs, z(:,ii+ninit-1)', W);
  x_0(:,1) = x_0(:,2);
  x_0(:,2) = x_0(:,3);
  x_0(:,3) = x(1:3,ii);
  P_0(:,:,1) = P_0(:,:,2);
  P_0(:,:,2) = P_0(:,:,3);
  P_0(:,:,3) = P(1:3,1:3,ii);
  [x(4:6,ii), P(4:6,4:6,ii)] = traj_pos2vel(dt, x_0(:,1),  P_0(:,:,1),...
                                                x_0(:,2),  P_0(:,:,2),...
                                                x_0(:,3),  P_0(:,:,3));
end
filt_rng = ninit:(ninit+nfilt-1);
res_plot("WLS", t(filt_rng), x_true(:,filt_rng), x, P);

