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
% filtering of a dynamic scenario.  A simple drag free trajectory is created
% along with simulated range observations (using range trackers from previous
% examples) that are subject to only measurement noise.  No bias effects
% have been added.  The goal is to implement a simple filter (observation model
% and dynamic model) that outperforms the WLS method (observation model only).
% This will be the basis for additional simulations as more measurement bias
% and systematic error are added.
%
% 1) Batch Method (observation model only).  A WLS method is used to
%    process range values at each time step.  The current position estimate
%    and two prior estimates, along with the corresponding covariances, are
%    then used to estimate the velocity and velocity covariance associated
%    with the latest position estimate.
% 
% 2) A linearized extended version of the U-D filter is implemented.
%
% Kurt Motekew  2016/09/28
%               2016/11/26
%

close all;
clear;

%
% Global Variable
%
% Set drag coefficient to zero for this initial filtering example
%

global global_b;
global_b = 0;

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
srng = .001;                                % Obs standard deviation
vrng = srng*srng;                           % Obs variance
W = vrng^-1;                                % Obs weighting "matrix"

t0 = 0;                                     % Scenario start time
tf = 3;                                     % Scenario stop time
dt = .01;                                   % Data rate

  % 'True' trajectory for which to generate simulated observations
[t, x_true] = traj_integ(t0, dt, tf, x0);
  % number of measurement sets - total obs = nsets*ntkrs
nsets = size(t,2);

%
% Create simulated range measurements
%

  % Vectors of measurement sets.  Columns at different times
z = zeros(ntkrs,nsets);
for ii = 1:nsets
  for jj = 1:ntkrs
    svec = x_true(1:3,ii) - tkrs(:,jj);
    z(jj,ii) = norm(svec) + srng*randn;
  end
end

  % Prime filter with batch derived state.  Use first three measurement
  % sets to generate three positions.  Then estimate velocity using the
  % three position estimates.  This creates an initial state for the third
  % time at which data is recorded.
ninit = 3;                                  % Number of initial batch estimates
x_0 = zeros(3, ninit);                      % Position vectors
P_0(3,3,ninit) = 0;                         % Position covariances
for ii = 1:ninit
  [x_0(:,ii), P_0(:,:,ii), ~] = box_locate(tkrs, z(:,ii)', W);
end
[vel, SigmaV] = traj_pos2vel(dt, x_0(:,1),  P_0(:,:,1),...
                                 x_0(:,2),  P_0(:,:,2),...
                                 x_0(:,3),  P_0(:,:,3));

  % Number of filtered outputs will be the initial batch derived
  % estimate plus the remaining number of observations
nfilt = nsets - ninit + 1;
x = zeros(6,nfilt);
P = zeros(6,6,nfilt);
x(:,1) = [x_0(:,ninit) ; vel]; 
P(1:3,1:3,1) = P_0(:,:,3);
P(4:6,4:6,1) = SigmaV;
  % Subset of full obs for filtering
filt_ndxoff = ninit - 1;
filt_rng = ninit:(filt_ndxoff + nfilt);

  %
  % Batch method - observation model only
  % Note offset in observation indexing
  %

for ii = 2:nfilt
  [x(1:3,ii), P(1:3,1:3,ii)] = box_locate(tkrs, z(:,ii+filt_ndxoff)', W);
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
res_plot('WLS', t(filt_rng), x_true(:,filt_rng), x, P);
  % Plot geometry
traj_plot(x, x_true(:,filt_rng), tkrs, blen);
title('WLS Trajectory');
view([70 20]);

  %
  % Linearized extended filter (via U-D method)
  % Nonlinear yet simple analytic trajectory for state propagation.
  % State vector of position and velocity - no acceleration terms in
  % the state transition matrix
  %

x_hat = x(:,1);                             % 'a priori' estimate and
P_hat = P(:,:,1);                           % covariance
[U, D] = mth_udut2(P_hat);                  % Decompose covariance for U-D form
x = zeros(6,nfilt);                         % Reset stored estimates
P = zeros(6,6,nfilt);
x(:,1) = x_hat;                             % Set first estimate to
P(:,:,1) = P_hat;                           % 'a priori' values
Phi = traj_strans(dt);                      % Fixed state transition matrix
for ii = 2:nfilt
    % First propagate state and covariance to new time - no process noise
  [~, U, D] = est_pred_ud(x_hat, U, D, Phi, zeros(6), eye(6));
  pos = traj_pos(dt, x_hat(1:3), x_hat(4:6));
  vel = traj_vel(dt, x_hat(4:6));
  x_bar = [pos ; vel];
  
    % Obs update based on observed (z) vs. computed (zc) residual (r)
  for jj = 1:ntkrs
    Ax = zeros(1,6);
    Ax(1:3) = est_drng_dloc(tkrs(:,jj), x(1:3,ii));
    zc = norm(x_bar(1:3) - tkrs(:,jj));
    r = z(jj,ii+filt_ndxoff) - zc;
    [x_hat, U, D] = est_upd_ud(x_bar, U, D, Ax, r, vrng);
    x_bar = x_hat;
  end
  P_hat = U*D*U';
  x(:,ii) = x_hat;
  P(:,:,ii) = P_hat;
end
res_plot('Linearized Extended UD', t(filt_rng), x_true(:,filt_rng), x, P);
  % Plot geometry
traj_plot(x, x_true(:,filt_rng), tkrs, blen);
title('Linearized Extended UD Trajectory');
view([70 20]);
