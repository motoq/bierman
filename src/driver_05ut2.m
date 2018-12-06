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
% Schmidt consider Kalman, SRIF, and UKF with consider parameter
% incorporation methods with process noise (prediction update)
% and bias (observation update error).  An error is added
% to the positions of each tracker.  Different random errors are added
% to each component of each tracker.  This bias is held constant for each
% trajectory.  This script adds the UKF so it can be compared to the SKF
% and SRIF methods.
%
% Note that unlike most of the other SKF examples in this project,
% the estimate/consider cross covariance is not scaled ("tuned").
%
% Also note the observation bias has been bumped up by a factor of 10!
%
% Kurt Motekew  2018/12/01
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
sblen = .01*blen;                          % .1% of total distance for tracker
SigmaBlen = sblen*sblen;
ny = 3;
tkrs = [                                    % uncertainty
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

  % Bias tracker locations along each axis - constant for each trajectory
tkr_deltas = sblen*randn(3,ntkrs);
  % Vectors of measurement sets.  Columns at different times
z = zeros(ntkrs,nsets);
for ii = 1:nsets
  for jj = 1:ntkrs
    svec = x_true(1:3,ii) - (tkrs(:,jj) + tkr_deltas(:,jj));
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
  % Linearized extended Schmidt Kalman filter
  %
x_hat = x_hat0;
P_hat = P_hat0;
x = zeros(6,nfilt);                         % Reset stored estimates
P = zeros(6,6,nfilt);
x(:,1) = x_hat;                             % Set first estimate to
P(:,:,1) = P_hat;                           % 'a priori' values
Phi = traj_strans(dt);                      % Fixed state transition matrix
SigmaZ = vrng;
SigmaY = SigmaBlen*ones(ny);
SigmaXY = zeros(6,ny);
tic;
for ii = 2:nfilt
    % First propagate state and covariance to new time - add
    % process noise when propagating covariance
  pos = traj_pos(dt, x_hat(1:3), x_hat(4:6));
  vel = traj_vel(dt, x_hat(4:6));
  x_bar = [pos ; vel];
  %Q = (.5*global_b)^2;
  Q = (2*global_b)^2;
  G = -[.5*vel*dt ; vel]*dt;
  SigmaX = Phi*P_hat*Phi' + G*Q*G';
    % Obs update based on observed (z) vs. computed (zc) residual (r)
  for jj = 1:ntkrs
    Ax = zeros(1,6);
    Ax(1:3) = est_drng_dloc(tkrs(:,jj), x_bar(1:3));
    Ay = est_drng_dpos(tkrs(:,jj), x_bar(1:3));
    zc = norm(x_bar(1:3) - tkrs(:,jj));
    r = z(jj,ii+filt_ndxoff) - zc;
    [x_bar, SigmaX, SigmaXY] = est_upd_schmidt(x_bar, SigmaX, Ax,...
                                                      SigmaY, Ay,...
                                                      SigmaXY, r, SigmaZ);
  end
  x_hat = x_bar;
  P_hat = SigmaX;
  x(:,ii) = x_hat;
  P(:,:,ii) = P_hat;
end
sck_time = toc;
res_plot('Linearized Extended Schmidt Consider Kalman',...
         t(filt_rng), x_true(:,filt_rng), x, P);
  % Plot geometry
traj_plot(x, x_true(:,filt_rng), tkrs, blen);
title('Linearized Extended Kalman Trajectory with Q & Y');
view([70 20]);

  %
  % SRIF
  %
x_hat = x_hat0;
P_hat = P_hat0;
ATA = P_hat^-1;
Rx = mth_sqrtm(ATA);
b = Rx*x_hat;                               % Lazy, properly size b of zeros
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
Rxy = zeros(6,ny);
Ay = zeros(ntkrs,ny);
Py0 = SigmaBlen*eye(ny);
Ry = sqrtm(Py0^-1);
by = zeros(ny,1);
r = zeros(ntkrs,1); 
tic;
for ii = 2:nfilt
    % First propagate state and information array to new time - add
    % process noise when propagating information array
  pos = traj_pos(dt, x_hat(1:3), x_hat(4:6));         
  vel = traj_vel(dt, x_hat(4:6));
  x_bar = [pos ; vel];           
  G = -[.5*vel*dt ; vel]*dt;
  [Rx, b] = est_pred_hhsrif(Rx, b, PhiInv, Rw, G);
  b = 0*b;
    % Obs update based on observed (z) vs. computed (zc) residual (r)
  for jj = 1:ntkrs
    Ax(jj,1:3) = est_drng_dloc(tkrs(:,jj), x_bar(1:3));
    Ay(jj,1:3) = est_drng_dpos(tkrs(:,jj), x_bar(1:3));
    zc = norm(x_bar(1:3) - tkrs(:,jj));
    r(jj) = z(jj,ii+filt_ndxoff) - zc;
  end
  [dx, Rx, Rxy, b, Ry, by] = est_upd_hhsrif_bias(Rx, Rxy, b, Ax, r, WsqrtSet,...
                                                                    Ry, by, Ay);
  b = 0*b;
  by = 0*by;
  x_hat = x_bar + dx;
    % Tapley, et. al. covariance method
  %Rn = [ Rx Rxy ; zeros(3,6) Ry ];
  %Rninv = mth_triinv(Rn);
  %Rxinv = Rninv(1:6,1:6);
  %S = -Rxinv*Rxy;
  %Pn = Rninv*Rninv';
  %Px = Pn(1:6,1:6);
  %Py = Pn(7:9,7:9);
  %P_hat = Px + S*Py*S';
    % Bierman covariance method
  Rxinv = mth_triinv(Rx);
  S = -Rxinv*Rxy;
  Phatc = Rxinv*Rxinv';
  Ryinv = mth_triinv(Ry);
  Py = Ryinv*Ryinv';
  P_hat = Phatc + S*Py*S';

  x(:,ii) = x_hat;
  P(:,:,ii) = P_hat;
end
srif_time = toc;
res_plot('Linearized Extended SRIF with Q & Y',...
         t(filt_rng), x_true(:,filt_rng), x, P);    
  % Plot geometry                               
traj_plot(x, x_true(:,filt_rng), tkrs, blen);
title('Linearized Extended SRIF Trajectory with Q & Y');
view([70 20]);


  %
  % Unscented Kalman with augmented state
  %

x_hat = x(:,1);                             % 'a priori' estimate and
P_hat = P(:,:,1);                           % covariance
x = zeros(6,nfilt);                         % Reset stored estimates
P = zeros(6,6,nfilt);
x(:,1) = x_hat;                             % Set first estimate to
P(:,:,1) = P_hat;                           % 'a priori' values
SigmaZ = srng*srng*eye(ntkrs);              % Observation covariance
alpha = .69;                                % Dude!
kappa = 0;
beta = 2;                                   % Gaussian

nx = size(x_hat,1);                         % # solve for params
nc = ny*ntkrs;                              % # total consider params
nx_a = nx + nc;                             % Augmented state vector size
SigmaXY = zeros(nx,nc);                     % Cross covariance
SigmaY = SigmaBlen*eye(nc);                 % Fixed consider covariance
x_hat_a = zeros(nx_a,1);                    % Augmented state vector
  % Weights for time and obs updates based on sigma vectors
tic;
for ii = 2:nfilt
    % Fill in augmented state vector and covariance
  x_hat_a(1:nx,1) = x_hat;
  for jj = 1:ntkrs
    ndx1 = nx + (jj-1)*ny + 1;
    ndx2 = ndx1 + 2;
    x_hat_a(ndx1:ndx2,1) = tkrs(:,jj);
  end
  P_hat_a = [P_hat SigmaXY ; SigmaXY' SigmaY];
    % Sigma vectors
  [Chi, w_m, w_c] = est_ut_sigma_vec(x_hat_a, P_hat_a, alpha, kappa, beta);
  dim = size(Chi,1);
  n_sigma_vec = size(Chi, 2);
    % Propagate sigma vectors - consider parameters remain the same
  for kk = 1:n_sigma_vec
    pos = traj_pos(dt, Chi(1:3,kk), Chi(4:6,kk));
    vel = traj_vel(dt, Chi(4:6,kk));
    Chi(1:3,kk) = pos;
    Chi(4:6,kk) = vel;
  end
    % Propagated estimate and covariance
  Q = (2*global_b)^2;
  G = -[.5*Chi(4:6,1)*dt ; Chi(4:6,1)]*dt;
  SigmaZq = G*Q*G';
  [x_bar_a, P_bar_a] = est_pred_ukf(Chi, w_m, w_c,...
                               [SigmaZq zeros(nx,nc) ; zeros(nc,nx) zeros(nc)]);    % Redraw sigma points to incorporate process noise effects
  [Chi, w_m, w_c] = est_ut_sigma_vec(x_bar_a, P_bar_a, alpha, kappa, beta);
    % Computed sigma vector based obs - note use of consider tracker locations
  Z = zeros(ntkrs,n_sigma_vec);
  for jj = 1:ntkrs
    ndx1 = nx + (jj-1)*ny + 1;
    ndx2 = ndx1 + 2;
    for kk = 1:n_sigma_vec
      Z(jj,kk) = norm(Chi(1:3,kk) - Chi(ndx1:ndx2,kk));
    end
  end
    % Update estimate based on available observations
  [x_hat_a, P_hat_a] = est_upd_ukf(x_bar_a, P_bar_a, Chi, w_m, w_c, Z,...
                                            z(:,ii+filt_ndxoff), SigmaZ);
    % Non-augmented state vector and covariance
    % Retain cross-covariance but forget consider parameter updates
  x_hat = x_hat_a(1:nx,1);
  P_hat = P_hat_a(1:nx,1:nx);
  SigmaXY = P_hat_a(1:nx,(nx+1):nx_a);

  x(:,ii) = x_hat;
  P(:,:,ii) = P_hat;
end
uskf_time = toc;
res_plot('USKF with Process Noise and Bias',...
          t(filt_rng), x_true(:,filt_rng), x, P);
  % Plot geometry
traj_plot(x, x_true(:,filt_rng), tkrs, blen);
title('USKF with Process Noise and Bias');
view([70 20]);

fprintf('\n SKF Time:\t%1.4f seconds', sck_time);
fprintf('\n SRIF Time:\t%1.4f seconds', srif_time);
fprintf('\n USKF Time:\t%1.4f seconds', uskf_time);
fprintf('\n');

