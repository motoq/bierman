%
% Copyright 2018 Kurt Motekew
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%

%
% Simple wiffle ball in a room trajectory
%
% This driver first estimats the trajectory using the batch method of
% the SRIF where all measurements taken at an instant in time are processed
% at once given range observations subject to noise while the filter's
% model is subject to process noise because it ignores drag (see driver_04b).
% Next, the true trajectory used as a reference trajectory while only the
% covariance is computed based on the error budget.
%
% Kurt Motekew  2018/03/21
%

close all;
clear;

%
% Global Variable
%

global global_b;
global_b = 0.05;

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
Wsqrt = 1/srng;
W = Wsqrt*Wsqrt';                           % Obs weighting "matrix"

t0 = 0;                                     % Scenario start time
tf = 3;                                     % Scenario stop time
dt = .01;                                   % Data rate

  % Create an ideal trajectory.  This is the model used by the filter
  % without the influence of observations.
[~, x_ideal] = traj_ideal(t0, dt, tf, x0);
  % 'True' trajectory for which to generate simulated observations
[t, x_true] = traj_integ(t0, dt, tf, x0);
  % number of measurement sets - total obs = nsets*ntkrs
nsets = size(t,2);

  % Plot geometry
traj_plot(x_ideal, x_true, tkrs, blen);
title('Boxed Tracker Geometry');
view([70 20]);

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
  % 'a priori' estimate and covariance to prime filters
  %

x_hat0 = x(:,1);                            % 'a priori' estimate and
P_hat0 = P(:,:,1);                          % covariance

  %
  % Householder Measurement Set Batch SRIF
  %

x_hat = x_hat0;
P_hat = P_hat0;
ATA = P_hat^-1;
R = mth_sqrtm(ATA);
b = R*x_hat;                                % Lazy, properly size b of zeros
b = 0*b;                                    % Differentials, not absolutes
x = zeros(6,nfilt);                         % Reset stored estimates
P = zeros(6,6,nfilt);
x(:,1) = x_hat;                             % Set first estimate to
P(:,:,1) = P_hat;                           % 'a priori' values
Phi = traj_strans(dt);                      % Fixed state transition matrix
PhiInv = Phi^-1;                            % Inverse for SRIF propagation
  % Fixed process noise
Q = diag((2*global_b)^2);
RwInv = sqrtm(Q);
Rw = RwInv^-1;
WsqrtSet = Wsqrt*eye(ntkrs);
Ax = zeros(ntkrs,6);
r = zeros(ntkrs,1);
for ii = 2:nfilt
    % First propagate state and information array to new time - add
    % process noise when propagating information array
  pos = traj_pos(dt, x_hat(1:3), x_hat(4:6));
  vel = traj_vel(dt, x_hat(4:6));
  x_bar = [pos ; vel];
  G = -[.5*vel*dt ; vel]*dt;
  [R, ~] = est_pred_hhsrif(R, b, PhiInv, Rw, G);
    % Obs update based on observed (z) vs. computed (zc) residual (r)
  for jj = 1:ntkrs
    Ax(jj,1:3) = est_drng_dloc(tkrs(:,jj), x_bar(1:3));
    zc = norm(x_bar(1:3) - tkrs(:,jj));
    r(jj) = z(jj,ii+filt_ndxoff) - zc;
  end
  [dx, R, ~, ~] = est_upd_hhsrif(R, b, Ax, r, WsqrtSet);
  x_hat = x_bar + dx;
  Rinv = mth_triinv(R);
  P_hat = Rinv*Rinv';
  x(:,ii) = x_hat;
  P(:,:,ii) = P_hat;
end
res_plot('Linearized Extended Batch SRIF with Q',...
         t(filt_rng), x_true(:,filt_rng), x, P);
  % Plot geometry
traj_plot(x, x_true(:,filt_rng), tkrs, blen);
title('Linearized Extended Batch SRIF Trajectory with Q');
view([70 20]);

P_hat = P_hat0;
ATA = P_hat^-1;
R = mth_sqrtm(ATA);
b = zeros(size(x_true(:,1), 1), 1);         % Differentials, not absolutes
P = zeros(6,6,nsets);
P(:,:,1) = P_hat;
Phi = traj_strans(dt);                      % Fixed state transition matrix
PhiInv = Phi^-1;                            % Inverse for SRIF propagation
  % Fixed process noise
Q = diag((2*global_b)^2);
RwInv = sqrtm(Q);
Rw = RwInv^-1;
WsqrtSet = Wsqrt*eye(ntkrs);
Ax = zeros(ntkrs,6);
r = zeros(ntkrs,1);
for ii = 2:nsets
    % Instead of propagating the state, pull the current reference trajectory
  vel = x_true(4:6,ii);
  G = -[.5*vel*dt ; vel]*dt;
  [R, ~] = est_pred_hhsrif(R, b, PhiInv, Rw, G);
    % Obs update based on observed (z) vs. computed (zc) residual (r)
  for jj = 1:ntkrs
    Ax(jj,1:3) = est_drng_dloc(tkrs(:,jj), x_true(1:3,ii));
  end
  [~, R, ~, ~] = est_upd_hhsrif(R, b, Ax, r, WsqrtSet);
  Rinv = mth_triinv(R);
  P_hat = Rinv*Rinv';
  P(:,:,ii) = P_hat;
end
cov_plot('Linearized Extended Batch SRIF with Q', t, P);
