clear;
close all;

%
% Illustrates linear least squares and constrained linear least squares
% fit to the curve y = a0 + a1*x + a2*exp(x) given an overdetermined
% system.
%
% The constrained fit example locks

%
% Define local functions.  For Matlab, move these three functions
% to the bottom of the file.  Octave requires them at the top of
% the file.
%

%
% Inputs:
%   x  independent variable
%   p  Coefficiens defining the curve, as defined in file header
%      p = [a0 a1 a2]
% Return:
%   function evaluated at x
%
function y = f(x, p)
  y = p(1) + p(2)*x + p(3)*exp(x);
end

%
% Inputs:
%   x  independent variable
%   p  Coefficiens defining the curve, as defined in file header
%      p = [a0 a1 a2]
% Return:
%   Derivative of function w.r.t. independent variable x
%
function yp =  dfdx(x, p)
  yp = p(2) + p(3)*exp(x);
end

%
% Input:
%   x  independent variable
% Return:
%   Derivative of function w.r.t. fit coefficients
%   [df/da0 df/da1 df/da2]
%
function Ap = dydp(x)
  Ap = [1 x exp(x)];
end

%
% Input:
%   x  independent variable
% Return:
%   Derivative of function w.r.t. fit coefficients w.r.t. independent
%   variables
%   [df2/da0dx df2/da1dx df2/da2dx]
%
function Axp = dydpdx(x)
  Axp = [0 1 exp(x)];
end

%
% End local functions
%

  % Solve-for parameters - curve fit parameters a0, a1, a2
np = 3;

  % yvals are the observations used for the curve fit
xvals = [0,  0.5, 1, 1.5, 2, 2.5,  3,  3.5,  4]';
yvals = [1,  2,   3, 5,   7, 8,   10, 14,    20]';
n = size(xvals,1);
  % LLS solution
Ap = zeros(n, np);
for ii = 1:n
  Ap(ii,:) = dydp(xvals(ii));
end
phat = ((Ap'*Ap)^-1)*Ap'*yvals;
  % Plot solution
xfit = 0:.2:4;
nfit = size(xfit,2);
yfit = zeros(1,nfit);
for ii = 1:nfit
  yfit(ii) = f(xfit(ii), phat);
end
figure;
plot(xfit, yfit, xvals, yvals, 'o')
xlabel('x');
ylabel('y');
title('LLS Fit');
%print -deps lls.eps

  % Constraints force endpoints to exactly match
  % and force a slope of 2 at the first point
slope = 0;
xcvals = [xvals(1), xvals(n), xvals(1)]';
ycvals = [yvals(1), yvals(n), slope]';
nc = size(xcvals, 1);

Ac = zeros(nc, np);
  % First set partials for endpoint constraints
for ii = 1:(nc-1)
  Ac(ii,:) = dydp(xcvals(ii));
end
  % Slope constraint
Ac(nc,:) = dydpdx(xcvals(nc));

  % Augmented normal equations
A = [2*Ap'*Ap  Ac' ; Ac zeros(nc, nc)];
yc = [2*Ap'*yvals ; ycvals];
  % Lagrange multipliers are the last three elements of the
  % augmented solve-for parameter vector
phat = ((A'*A)^-1)*A'*yc;

for ii = 1:nfit
  yfit(ii) = f(xfit(ii), phat);
end
figure;
plot(xfit, yfit, xvals, yvals, 'o');
xlabel('x');
ylabel('y');
title('Constrained LLS Fit');

  % Verify constraints are met
fprintf('\ny = %1.3f vs. y = %1.3f constraint at x = %1.1f',...
        f(xcvals(1), phat), ycvals(1), xcvals(1));
fprintf('\ny = %1.3f vs. y = %1.3f constraint at x = %1.1f',...
        f(xcvals(2), phat), ycvals(2), xcvals(2));
fprintf('\ndydx = %1.3f vs. dydx = %1.3f constraint at x = %1.1f',...
        dfdx(xcvals(3), phat), slope, xcvals(3));

fprintf('\n');
