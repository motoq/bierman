% Illustrates sequential estimation via Householder SRIF for a linear problem
% involving fitting a number of data points to the curve y = a + b*x + c*e^x
% where [a b c]' are the solve for parameters.
%
% The "truth" fit is accomplished using Householder triangularization of the
% information array.  The fit of the coefficients are then plotted against the
% data points.
%
% Next, the minimum number of data points are used to form an initial batch
% estimate.
%
% Finally, the remaining observations are used to refine the estimate through
% the sequential form of the Householder triangularzation method.
%
% Kurt Motekew   2016/08/23
%


clear all;
close all;

  % Data
x = [0; 1; 2; 3; 4];
y = [1; 3; 7; 10; 20];
sigma = 0.5;

m = 5;
n = 3;

  % Batch via Householder
A = zeros(m, 3);
for ii = 1:m
  A(ii, 1) = 1;
  A(ii, 2) = x(ii);
  A(ii, 3) = exp(x(ii));
end
SqrtW = eye(m)/sigma;
  % Info Array and decomposition
Rz0 = [ SqrtW*A SqrtW*y ];
Rz = mth_householder_tri(Rz0, 1);
R = Rz(1:n,1:n);
z = Rz(1:n,(n+1));
  % Coefficients
phat = mth_trisol(R,z);

% Plot line
start = 0;
incre = .2;
stop  = 4;
xc = start:incre:stop;
xc = xc';
yc = xc;
num = (stop-start)/incre + 1;
for ii = 1:num
  yc(ii) = phat(1) + phat(2)*xc(ii) + phat(3)*exp(xc(ii));
end
figure; hold;
plot(x, y, 'o', xc, yc, '-b');
%axis([-0.1, 4.1, -.5, 20.5]);
xlabel('x');
ylabel('y');
title('Curve Fit via Householder Triangularization of the Information Array');
%print -deps lls.eps;


  % Batch init via Householder
m = 3;
m2 = 5;

  % Batch via Householder
A = zeros(m, 3);
for ii = 1:m
  A(ii, 1) = 1;
  A(ii, 2) = x(ii);
  A(ii, 3) = exp(x(ii));
end
SqrtW = eye(m)/sigma;
  % Info Array and decomposition
Rz0 = [ SqrtW*A SqrtW*y(1:m,1) ];
Rz = mth_householder_tri(Rz0, 1);
R = Rz(1:n,1:n);
z = Rz(1:n,(n+1));
  % Coefficients
phat = mth_trisol(R,z);
  % Plot initial estimate
for ii = 1:num
  yc(ii) = phat(1) + phat(2)*xc(ii) + phat(3)*exp(xc(ii));
end
plot(xc, yc, '-r');

sw = 1/sigma;
for ii = (m+1):m2
  Rz0 = [R z ; sw*[1 x(ii) exp(x(ii))] sw*y(ii,1) ];
  Rz = mth_householder_tri(Rz0, 1);
  R = Rz(1:n,1:n);
  z = Rz(1:n,(n+1));
  phat = mth_trisol(R,z);
  % Need to fix 'legend' command if uncommenting
  %for jj = 1:num
  %  yc(jj) = phat(1) + phat(2)*xc(jj) + phat(3)*exp(xc(jj));
  %end
  %plot(xc, yc, '.r');
end
  % Plot updated estimate
for ii = 1:num
  yc(ii) = phat(1) + phat(2)*xc(ii) + phat(3)*exp(xc(ii));
end
plot(xc, yc, '*g');

legend('Data', 'Batch', 'Init', 'Seq');
