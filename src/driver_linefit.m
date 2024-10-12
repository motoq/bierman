%
% This driver script does a line fit to noisy data using standard WLS
% via normal equations and compares it to "squaring up" of the observation
% model via Householder transformations.  Separate plots are made for each
% method.  The normal equations method also plots 3-sigma bounds based on the
% estimate covariance.  Estimate and covariance differences between the
% normal equations method and Householder are printed to stdout.
%
% The full decomposition approach to the Householder reflection method
% is also illustrated using backwards substitution for the solve-for
% parameters and inversion of the upper triangular square root
% information matrix.
%
% Kurt Motekew 2024/10/11
% 

close all;
clear;

  % Number of standard deviations for which to plot covariance bounds
nsigma = 3;

  % Truth Values: y = ax + b
a = 0.3;                          % Slope
b = 1.0;                          % y-intercept

  % Y values are measurements - random uncertainty to apply - 1 sigma Gaussian
sigma_y = 0.05;

  % Generate simulated y values over given x values:  y = ax + b
x = (1:.01:2)';
m = size(x,1);
y = a*x + b + sigma_y*randn(m,1);

  % Jacobian - partials of obs w.r.t. solve-for
A = zeros(m,2);
A(:,1) = x;
A(:,2) = 1;

  % Observation covariance
Cy = sigma_y*sigma_y*eye(m);
  % Weighting matrix
W = Cy^-1;

  % Form normal equations for standard WLS method
ATWA = A'*W*A;
atwy = A'*W*y;
  % Estimate covariance via inverse of information matrix
Cp = ATWA^-1;
  % Estimate - just use information  matrix inverse vs.
  % something more numerically robust for this example
phat = Cp*atwy;

  % Plot line using estimated parameters...
xhat = 0:.1:2;
yhat = phat(1)*xhat + phat(2);
  % and uncertainty bounds
sigma_a = sqrt(Cp(1,1));
sigma_b = sqrt(Cp(2,2));
da = nsigma*sigma_a;
db = nsigma*sigma_b;
  % Jut plot upper limint, lower, increased and decreased slopes
yhatu = phat(1)*xhat + phat(2)+db;
yhatl = phat(1)*xhat + phat(2)-db;
yhati = (phat(1)+da)*xhat + phat(2);
yhatd = (phat(1)-da)*xhat + phat(2);

figure;  hold on;
plot(x, y, '.k');
plot(xhat, yhat, '-b');
plot(xhat, yhatu, '-r');
plot(xhat, yhatl, '-r');
plot(xhat, yhati, '-m');
plot(xhat, yhatd, '-m');
xlabel('x');
ylabel('y');
title('Normal Equations Solution');
axis equal;

  % Solve using Householder transformation of observation model

  % Normalize observation covariance such that Cy = I
S = chol(W);
  % Condense information, squaring up observation model
Ay = [S*A S*y];
R = mth_householder_tri(Ay, 1);
  % Backwards substitution for estimate and form covariance
Rp = R(1:2, 1:2);
phat_h = mth_trisol(Rp, R(1:2,3));
Sp = mth_triinv(Rp);
Cp_h = Sp*Sp';

  % Plot line using estimated parameters...
yhat = phat_h(1)*xhat + phat_h(2);
figure;  hold on;
plot(x, y, '.k');
plot(xhat, yhat, '-b');
xlabel('x');
ylabel('y');
title('Square Root Information Matrix Solution');
axis equal;

  % Summary info
fprintf('\nDifference between estimates: %1.3e', norm(phat - phat_h));
fprintf('\nDifference between ATWA and SA covariance matrices: %1.3e',...
        norm(norm(Cp - Cp_h)));

fprintf('\n');



 
