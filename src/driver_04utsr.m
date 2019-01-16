%
% Copyright 2019 Kurt Motekew
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%

%
% Simple wiffle ball in a room trajectory
%
% This driver is yet another extension of the _04 series comparing
% the extended Kalman to the unscented Kalman and a square root version of
% the unscented Kalman where both random observation noise and process
% "noise" is present.
%
% Kurt Motekew  2019/01/11
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
  % Extended filter via Kalman stabilized method with covariance
  % inflation.
  %

x_hat = x_hat0;
P_hat = P_hat0;
x = zeros(6,nfilt);                         % Reset stored estimates
P = zeros(6,6,nfilt);
x(:,1) = x_hat;                             % Set first estimate to
P(:,:,1) = P_hat;                           % 'a priori' values
Phi = traj_strans(dt);                      % Fixed state transition matrix
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
  P_bar = Phi*P_hat*Phi' + G*Q*G';
    % Obs update based on observed (z) vs. computed (zc) residual (r)
  for jj = 1:ntkrs
    Ax = zeros(1,6);
    Ax(1:3) = est_drng_dloc(tkrs(:,jj), x_bar(1:3));
    zc = norm(x_bar(1:3) - tkrs(:,jj));
    r = z(jj,ii+filt_ndxoff) - zc;
    [x_bar, P_bar] = est_upd_kalman(x_bar, P_bar, Ax, r, Wsqrt);
  end
  x_hat = x_bar;
  P_hat = P_bar;
  x(:,ii) = x_hat;
  P(:,:,ii) = P_hat;
end
kal_time = toc;
res_plot('EKF with Process Noise',...
         t(filt_rng), x_true(:,filt_rng), x, P);
  % Plot geometry
traj_plot(x, x_true(:,filt_rng), tkrs, blen);
title('EKF with Process Noise');
view([70 20]);

  %
  % Unscented Kalman version.
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
  % Weights for time and obs updates based on sigma vectors
tic;
for ii = 2:nfilt
    % Sigma vectors
  [Chi, w_m, w_c] = est_ut_sigma_vec(x_hat, P_hat, alpha, kappa, beta);
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
  Q = (2*global_b)^2;
  G = -[.5*Chi(4:6,1)*dt ; Chi(4:6,1)]*dt;
  SigmaZq = G*Q*G';
  [x_bar, P_bar] = est_pred_ukf(Chi, w_m, w_c, SigmaZq);
    % Redraw sigma points to incorporate process noise effects
  [Chi, w_m, w_c] = est_ut_sigma_vec(x_bar, P_bar, alpha, kappa, beta);
    % Computed sigma vector based obs                   
  Z = zeros(ntkrs,n_sigma_vec);
  for jj = 1:ntkrs
    for kk = 1:n_sigma_vec
      Z(jj,kk) = norm(Chi(1:3,kk) - tkrs(:,jj));
    end
  end  
    % Update estimate based on available observations
  [x_hat, P_hat] = est_upd_ukf(x_bar, P_bar, Chi, w_m, w_c, Z,...
                                             z(:,ii+filt_ndxoff), SigmaZ);
  x(:,ii) = x_hat;                                                        
  P(:,:,ii) = P_hat;
end
ukf_time = toc;
res_plot('UKF with Process Noise',...
         t(filt_rng), x_true(:,filt_rng), x, P);
  % Plot geometry
traj_plot(x, x_true(:,filt_rng), tkrs, blen);
title('UKF with Process Noise');
view([70 20]);

  %
  % Square Root Unscented Kalman filter
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
  % Weights for time and obs updates based on sigma vectors
[w_m, w_c] = est_ut_sigma_weights(size(x_hat,1), alpha, kappa, beta);
if w_c(1) < 0
  fprintf('\nNegative 0th Covariance Weighting not Supported\n');
  return;
end
sr_w_c = sqrt(w_c);
S_hat = mth_sqrtm(P_hat);
tic;
for ii = 2:nfilt
  Chi = est_ut_srsigma_vec(x_hat, S_hat, alpha, kappa);
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
  [x_bar, L_bar] = est_pred_srukf(Chi, w_m, sr_w_c, G*srQ);
    % Redraw sigma points to incorporate process noise effects
  Chi = est_ut_srsigma_vec(x_bar, L_bar, alpha, kappa);
    % Computed sigma vector based obs
  Z = zeros(ntkrs,n_sigma_vec);
  for jj = 1:ntkrs
    for kk = 1:n_sigma_vec
      Z(jj,kk) = norm(Chi(1:3,kk) - tkrs(:,jj));
    end
  end
    % Update estimate based on available observations
  [x_hat, S_hat] = est_upd_srukf(x_bar, L_bar, Chi, w_m, sr_w_c, Z,...
                                               z(:,ii+filt_ndxoff), SrZ);
  P_hat = S_hat*S_hat';
  
  x(:,ii) = x_hat;
  P(:,:,ii) = P_hat;
end
srukf_time = toc;
res_plot('SRUKF', t(filt_rng), x_true(:,filt_rng), x, P);
  % Plot geometry
traj_plot(x, x_true(:,filt_rng), tkrs, blen);
title('SRUKF Trajectory');
view([70 20]);

fprintf('\n Kalman Time:\t%1.4f seconds', kal_time);
fprintf('\n UKF Time:\t\t%1.4f seconds', ukf_time);
fprintf('\n SRUKF Time:\t\t%1.4f seconds', srukf_time);
fprintf('\n');

