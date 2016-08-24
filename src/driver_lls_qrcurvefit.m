% Illustrates sequential estimation via QR SRIF for a linear problem
% involving fitting a number of data points to the curve y = a + b*x + c*e^x
% where [a b c]' are the solve for parameters.
%
% The "truth" fit is accomplished using QR decomposition of the
% Jacobian.  The fit of the coefficients are then plotted against the
% data points.
%
% Next, the minimum number of data points are used to form an initial batch
% estimate.
%
% Finally, the remaining observations are used to refine the estimate through
% the sequential form of the QR SRIF method.
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

  % Batch via QR
A = zeros(m, 3);
for ii = 1:m
  A(ii, 1) = 1;
  A(ii, 2) = x(ii);
  A(ii, 3) = exp(x(ii));
end
SqrtW = eye(m)/sigma;
A = SqrtW*A;
z = SqrtW*y;
[Q, R] = mth_qr(A);
qty = Q'*z;
  % Coefficients
phat = mth_trisol(R, qty);

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
title('Curve Fit via QR decomposition of the Jacobian');
%print -deps lls.eps;


  % Batch init via QR
m = 3;
m2 = 5;

  % Batch via QR
A = zeros(m, 3);
for ii = 1:m
  A(ii, 1) = 1;
  A(ii, 2) = x(ii);
  A(ii, 3) = exp(x(ii));
end
SqrtW = eye(m)/sigma;
A = SqrtW*A;
z = SqrtW*y(1:m,1);
[Q, R] = mth_qr(A);
qty = Q'*z;
  % Coefficients
phat = mth_trisol(R, qty);
  % Plot initial estimate
for ii = 1:num
  yc(ii) = phat(1) + phat(2)*xc(ii) + phat(3)*exp(xc(ii));
end
plot(xc, yc, '-r');

sw = 1/sigma;
for ii = (m+1):m2
  A = [R ; sw*[1 x(ii) exp(x(ii))]];
  z = [qty ; sw*y(ii,1)];
  [Q, R] = mth_qr(A);
  qty = Q'*z;
  phat = mth_trisol(R, qty);
end
  % Plot updated estimate
for ii = 1:num
  yc(ii) = phat(1) + phat(2)*xc(ii) + phat(3)*exp(xc(ii));
end
plot(xc, yc, '*g');

legend('Data', 'Batch', 'Init', 'Seq');
