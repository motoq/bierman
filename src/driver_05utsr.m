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
% Schmidt Kalman and unscented Kalman with an augmented state vector.
% Process noise (prediction update) and bias (observation update error)
% are included and accounted for by each estimator.  An error is added
% to the positions of each tracker.  Different random errors are added
% to each component of each tracker.  This bias is held constant for each
% trajectory.  The UKF version accounts for bias by including the consider
% parameters in the augmented state vector while not actually updating
% them or their associated covariance.  For some reason the literature
% describes this as an unscented Schmidt Kalman filter, but it in no
% way resembles one.
%
% Note that unlike the other SKF examples in this project, the
% estimate/consider cross covariance is not scaled ("tuned").
%
% Square root version in progress...
%
% Kurt Motekew  2018/11/24
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
sblen = .001*blen;                          % .1% of total distance for tracker
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
  %SigmaXY = 0.7*SigmaXY;
  x(:,ii) = x_hat;
  P(:,:,ii) = P_hat;
end
sck_time = toc;
res_plot('SKF with Process Noise and Bias',...
         t(filt_rng), x_true(:,filt_rng), x, P);
  % Plot geometry
traj_plot(x, x_true(:,filt_rng), tkrs, blen);
title('SKF with Process Noise and Bias');
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
                               [SigmaZq zeros(nx,nc) ; zeros(nc,nx) zeros(nc)]);
    % Redraw sigma points to incorporate process noise effects
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

  %
  % Square Root Unscented Kalman filter with augmented state
  %

x_hat = x(:,1);                             % 'a priori' estimate and
P_hat = P(:,:,1);                           % covariance
x = zeros(6,nfilt);                         % Reset stored estimates
P = zeros(6,6,nfilt);
x(:,1) = x_hat;                             % Set first estimate to
P(:,:,1) = P_hat;                           % 'a priori' values
SigmaZ = srng*srng*eye(ntkrs);              % Observation covariance
SrZ = mth_sqrtm(SigmaZ);
alpha = .69;
kappa = 0;
beta = 2;                                   % Gaussian

nx = size(x_hat,1);                         % # solve for params
nc = ny*ntkrs;                              % # total consider params
nx_a = nx + nc;                             % Augmented state vector size
SigmaXY = zeros(nx,nc);                     % Cross covariance
SigmaY = SigmaBlen*eye(nc);                 % Fixed consider covariance
x_hat_a = zeros(nx_a,1);                    % Augmented state vector
x_hat_a(1:nx) = x_hat;
for ii = 1:ntkrs
  ndx1 = nx + (ii-1)*ny + 1;
  ndx2 = ndx1 + 2;
  x_hat_a(ndx1:ndx2,1) = tkrs(:,ii);
end
P_hat_a = [P_hat SigmaXY ; SigmaXY' SigmaY];
  % Get weights
[w_m, w_c] = est_ut_sigma_weights(nx_a, alpha, kappa, beta);
if w_c(1) < 0
  fprintf('\nNegative 0th Covariance Weighting not Supported\n');
  return;
end
sr_w_c = sqrt(w_c);
S_hat_a = mth_sqrtm(P_hat_a);
tic;
for ii = 2:nfilt
  Chi = est_ut_srsigma_vec(x_hat_a, S_hat_a, alpha, kappa);
  dim = size(Chi,1);
  n_sigma_vec = size(Chi, 2);
    % Propagate sigma vectors
  for kk = 1:n_sigma_vec
    pos = traj_pos(dt, Chi(1:3,kk), Chi(4:6,kk));
    vel = traj_vel(dt, Chi(4:6,kk));
    Chi(1:3,kk) = pos;
    Chi(4:6,kk) = vel;
  end
    % Propagated estimate and covariance
  srQ = 2*global_b;
  G = -[.5*Chi(4:6,1)*dt ; Chi(4:6,1)]*dt;
  %[x_bar_a, L_bar_a] = est_pred_srukf(Chi, w_m, sr_w_c, G*srQ);
  [x_bar_a, L_bar_a] = est_pred_srukf(Chi, w_m, sr_w_c, [G*srQ ; zeros(nc,1)]);
    % Redraw sigma points to incorporate process noise effects
  Chi = est_ut_srsigma_vec(x_bar_a, L_bar_a, alpha, kappa);
  
    % Computed sigma vector based obs
  Z = zeros(ntkrs,n_sigma_vec);
  for jj = 1:ntkrs
    ndx1 = nx + (jj-1)*ny + 1;
    ndx2 = ndx1 + 2;
    for kk = 1:n_sigma_vec
      Z(jj,kk) = norm(Chi(1:3,kk) - Chi(ndx1:ndx2,kk));
    end
  end
    % Update estimate based on available observations
  [x_hat_a, S_hat_a] = est_upd_srukf(x_bar_a, L_bar_a, Chi, w_m, sr_w_c, Z,...
                                              z(:,ii+filt_ndxoff), SrZ, nc);
  S_hat = S_hat_a(1:nx,1:nx);
  P_hat = S_hat*S_hat';
  
  x(:,ii) = x_hat_a(1:nx);
  P(:,:,ii) = P_hat;
end
srukf_time = toc;
res_plot('SRUKF', t(filt_rng), x_true(:,filt_rng), x, P);
  % Plot geometry
traj_plot(x, x_true(:,filt_rng), tkrs, blen);
title('SRUKF Trajectory');
view([70 20]);

fprintf('\n Schmidt EKF Time:\t%1.4f seconds', sck_time);
fprintf('\n Schmidt UKF Time:\t%1.4f seconds', uskf_time);
fprintf('\n SRUKF Time:\t\t%1.4f seconds', srukf_time);
fprintf('\n');

