%
% Copyright 2018 Kurt Motekew
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%

%
% This driver script illustrates propagation of a 3D covariance via the
% unscented transformation.  It is analogis to a position vector and
% covariance propagated in time where the velocity uncertainty and
% propagation model are perfect.
%
% A helix is used to represent the true trajectory.  Covariance analysis
% via the UT is performed.
%

clear;
close all;

  % Helix definition and range over which to plot
r = 10;
c = 1;
t_curve = 0:.1:(2*pi);

  % Plot Helix - truth trajectory
[x_curve, y_curve, z_curve] = mth_helix_parametric(t_curve, r, c);
figure; hold on;
plot3(x_curve, y_curve, z_curve);

  % Position estimate (on reference trajectory) and covariance
phat0 = [ x_curve(1) y_curve(1) z_curve(1)]';
SigmaP0 = [ 3 2 1 ; 2 4 3 ; 1 3 5 ];

  % Sigma Vectors at start time
alpha = .69;                           % Small positive value, 1e-4 <= a <= 1
kappa = 0; %3 - L;                     % Secondary scaling parameter
beta  = 2;
[Chi, w_m, w_c] = est_ut_sigma_vec(phat0, SigmaP0, alpha, kappa, beta);
n_params = size(Chi, 1);
n_sigma_vec = size(Chi, 2);
  % Plot sigma vectors and initial covariance
[XX, YY, ZZ] = matrix3X3_points(SigmaP0, 20);
mesh(XX + x_curve(1), YY + y_curve(1), ZZ + z_curve(1));
colormap([.1 .2 .3]);
hidden('off');
ChiPlt = Chi';
scatter3(ChiPlt(:,1), ChiPlt(:,2), ChiPlt(:,3), 'b', 'filled');

  % Compute 5 points for 4th order poly fit to helix
t_curve = 0:1:4;
[x_curve, y_curve, z_curve] = mth_helix_parametric(t_curve, r, c);
scatter3(x_curve, y_curve, z_curve);
ax = mth_lpoly_fit(t_curve, x_curve, 4);
ay = mth_lpoly_fit(t_curve, y_curve, 4);
az = mth_lpoly_fit(t_curve, z_curve, 4);
  % Propagate to 3rd (middle) point and plot
ii = size(t_curve,2);
ii = floor(ii/2 + 1);
t_prop = t_curve(ii);
x_cov = mth_lpoly_eval(ax, t_prop);
y_cov = mth_lpoly_eval(ay, t_prop);
z_cov = mth_lpoly_eval(az, t_prop);
scatter3(x_cov, y_cov, z_cov, 'filled');
  % Propagate sigma vectors
  % Modify "intercept" portion of poly def to match sigma vectors
axi = ax;
ayi = ay;
azi = az;
t = 0:.1:t_prop;
for ii = 1:n_sigma_vec
  axi(1) = Chi(1,ii);
  ayi(1) = Chi(2,ii);
  azi(1) = Chi(3,ii);
  Chi(1,ii) =  mth_lpoly_eval(axi, t_prop);
  Chi(2,ii) =  mth_lpoly_eval(ayi, t_prop);
  Chi(3,ii) =  mth_lpoly_eval(azi, t_prop);
  x_curve = mth_lpoly_eval(axi, t);
  y_curve = mth_lpoly_eval(ayi, t);
  z_curve = mth_lpoly_eval(azi, t);
  plot3(x_curve, y_curve, z_curve, 'b');
end
  % Plot propagated covariance at new location
[phat, SigmaP] = est_pred_ukf(Chi, w_m, w_c, zeros(n_params));
[XX, YY, ZZ] = matrix3X3_points(SigmaP, 20);
mesh(XX + phat(1), YY + phat(2), ZZ + phat(3));
hidden('off');
ChiPlt = Chi';
scatter3(ChiPlt(:,1), ChiPlt(:,2), ChiPlt(:,3), 'b', 'filled');

  % Plot polynomial approximation to helix
t_curve = 0:.1:4;
x_curve = mth_lpoly_eval(ax, t_curve);
y_curve = mth_lpoly_eval(ay, t_curve);
z_curve = mth_lpoly_eval(az, t_curve);
plot3(x_curve, y_curve, z_curve, 'r', 'LineWidth', 3)

xlabel('x');
ylabel('y');
zlabel('z');
title('Propagation of Position Covariance with No Velocity Uncertainty');
axis equal;
