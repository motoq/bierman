%
% Copyright 2018 Kurt Motekew
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%

%
% This driver script illustrates the unscented transformation of a 3D estimate
% given an associated covariance.
%

clear;
close all;

  % Estimate (origin of reference frame) and covariance
phat0 = [ 0 0 0]';
SigmaP0 = [ 3 2 1 ; 2 4 3 ; 1 3 5 ];

alpha = .69;                           % Small positive value, 1e-4 <= a <= 1
kappa = 0; %3 - L;                     % Secondary scaling parameter
beta  = 2;

[Chi, w_m, w_c] = est_ut_sigma_vec(phat0, SigmaP0, alpha, kappa, beta);
n_params = size(Chi, 1)
n_sigma_vec = size(Chi, 2);

  % Plot covariance and sigma vectors
figure; hold on;
matrix3X3_plot(SigmaP0, 20, false);
ChiPlt = Chi';
scatter3(ChiPlt(:,1), ChiPlt(:,2), ChiPlt(:,3), 'r', 'filled');
title('Unscented Transformation Sigma Vectors');




  % Helix definition and range over which to plot
r = 10;
c = 1;
t_curve = 0:.1:(2*pi);

[x_curve, y_curve, z_curve] = mth_helix_parametric(t_curve, r, c);
figure; hold on;
plot3(x_curve, y_curve, z_curve);

[XX, YY, ZZ] = matrix3X3_points(SigmaP0, 20);
mesh(XX + x_curve(1), YY + y_curve(1), ZZ + z_curve(1));
colormap([.1 .2 .3]);
hidden('off');

t_curve = 0:.3:1.2;
[x_curve, y_curve, z_curve] = mth_helix_parametric(t_curve, r, c);
scatter3(x_curve, y_curve, z_curve);
ax = mth_lpoly_fit(t_curve, x_curve, 4);
ay = mth_lpoly_fit(t_curve, y_curve, 4);
az = mth_lpoly_fit(t_curve, z_curve, 4);

ii = size(t_curve,2);
ii = floor(ii/2 + 1);
t_cov = t_curve(ii);
x_cov = mth_lpoly_eval(ax, t_cov);
y_cov = mth_lpoly_eval(ay, t_cov);
z_cov = mth_lpoly_eval(az, t_cov);
scatter3(x_cov, y_cov, z_cov, 'filled');

axi = ax;
ayi = ay;
azi = az;
for ii = 1:n_sigma_vec
  axi(1) = ax(1) + Chi(1,ii);
  ayi(1) = ay(1) + Chi(2,ii);
  azi(1) = az(1) + Chi(3,ii);
  Chi(1,ii) =  mth_lpoly_eval(axi, t_cov);
  Chi(2,ii) =  mth_lpoly_eval(ayi, t_cov);
  Chi(3,ii) =  mth_lpoly_eval(azi, t_cov);
end
[phat, SigmaP] = est_pred_ukf(Chi, w_m, w_c, zeros(n_params));
[XX, YY, ZZ] = matrix3X3_points(SigmaP, 20);
mesh(XX + phat(1), YY + phat(2), ZZ + phat(3));
hidden('off');

t_curve = 0:.1:1.2;
x_curve = mth_lpoly_eval(ax, t_curve);
y_curve = mth_lpoly_eval(ay, t_curve);
z_curve = mth_lpoly_eval(az, t_curve);
plot3(x_curve, y_curve, z_curve, 'r', 'LineWidth', 3)



xlabel('x');
ylabel('y');
zlabel('z');
axis equal;

%figure; hold on;
%matrix3X3_plot(SigmaP0, 20, false);
%matrix3X3_plot(SigmaP, 20, false);
