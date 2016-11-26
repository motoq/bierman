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
% 2) Adding first filtered versions - in progress
%
% Kurt Motekew  2016/09/28
%               2016/11/25
%

close all;
clear;

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
  % Subset of full obs for filtering
filt_rng = ninit:(ninit+nfilt-1);

  % Batch method - observation model only
  % Note offset in observations being passed
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
res_plot('WLS', t(filt_rng), x_true(:,filt_rng), x, P);
  % Plot geometry
traj_plot(x, x_true(:,filt_rng), tkrs, blen);
title('WLS Trajectory');
view([70 20]);

  % Linear filter (via U-D method)
Phi = traj_strans(dt);
x_hat = x(:,1);
x_bar = x_hat;
P_hat = P(:,:,1);
x = zeros(6,nfilt);
P = zeros(6,6,nfilt);
x(:,1) = x_hat;
P(:,:,1) = P_hat;
for ii = 2:nfilt
    % First propagate to new time
  [U, D] = mth_udut2(P_hat);
  [x_bar, U, D] = est_pred_ud(x_hat, U, D, Phi, zeros(6), eye(6));
  
    % Obs update
  for jj = 1:ntkrs
    Ap = zeros(1,6);
    Ap(1:3) = est_drng_dloc(tkrs(:,jj), x(1:3,ii));
    zc = norm(x_bar(1:3) - tkrs(:,jj));
    r = z(jj,ii+ninit-1) - zc;
    [x_bar, U, D] = est_upd_ud(x_bar, U, D, Ap, r, vrng);
  end
  x_hat = x_bar;
  P_hat = U*D*U';
  x(:,ii) = x_hat;
  P(:,:,ii) = P_hat;
end
res_plot('UD', t(filt_rng), x_true(:,filt_rng), x, P);
  % Plot geometry
traj_plot(x, x_true(:,filt_rng), tkrs, blen);
title('Linear UD Trajectory');
view([70 20]);

  % Linearized filter (via U-D method)
Phi = traj_strans(dt);
x_hat = x(:,1);
x_bar = x_hat;
P_hat = P(:,:,1);
x = zeros(6,nfilt);
P = zeros(6,6,nfilt);
x(:,1) = x_hat;
P(:,:,1) = P_hat;
for ii = 2:nfilt
    % First propagate to new time
  [U, D] = mth_udut2(P_hat);
  [~, U, D] = est_pred_ud(x_hat, U, D, Phi, zeros(6), eye(6));
  x_bar(1:3) = traj_pos(dt, x_hat(1:3), x_hat(4:6));
  x_bar(4:6) = traj_vel(dt, x_hat(4:6));
  
    % Obs update
  for jj = 1:ntkrs
    Ap = zeros(1,6);
    Ap(1:3) = est_drng_dloc(tkrs(:,jj), x(1:3,ii));
    zc = norm(x_bar(1:3) - tkrs(:,jj));
    r = z(jj,ii+ninit-1) - zc;
    [x_bar, U, D] = est_upd_ud(x_bar, U, D, Ap, r, vrng);
  end
  x_hat = x_bar;
  P_hat = U*D*U';
  x(:,ii) = x_hat;
  P(:,:,ii) = P_hat;
end
res_plot('Linearized UD', t(filt_rng), x_true(:,filt_rng), x, P);
  % Plot geometry
traj_plot(x, x_true(:,filt_rng), tkrs, blen);
title('Linear UD Trajectory');
view([70 20]);
