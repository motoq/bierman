function res_plot(title_hdr,t, x, xhat, SigmaX)
% RES_PLOT plots position and velocity error along with the associated 95%
% covariance.
%
%-----------------------------------------------------------------------
% Copyright 2016 Kurt Motekew
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%-----------------------------------------------------------------------
%
% Inputs
%   title_hdr   String indicating what method is used
%   t           Array of times
%   x           True position & velocity, [6xN]
%   xhat        Estimated position & velocity, [6xN]
%   SigmaX      xhat covariance, [6x6xN]
%
% Author:  Kurt Motekew    20160816
%

    % Scale covariance by 95% containment
  SF95_3D = 2.796;
  SigmaX = SF95_3D*SF95_3D*SigmaX;
  residual = abs(x - xhat);

    % Position
  figure; hold on;
  plot(t, residual(1,:), 'ro', t, residual(2,:), 'go', t, residual(3,:), 'bo');
  legend('x', 'y', 'z');
  plot(t, sqrt(squeeze(SigmaX(1,1,:))), 'r-',...
       t, sqrt(squeeze(SigmaX(2,2,:))), 'g-',...
       t, sqrt(squeeze(SigmaX(3,3,:))), 'b-');
  xlabel('Time');
  ylabel('Position Error, per Axis');
  title(strcat(title_hdr, ' Estimated Position Error and 95% Error Bounds'));

    % Velocity
  figure; hold on;
  plot(t, residual(4,:), 'ro', t, residual(5,:), 'go', t, residual(6,:), 'bo');
  legend('x', 'y', 'z');
  plot(t, sqrt(squeeze(SigmaX(4,4,:))), 'r-',...
       t, sqrt(squeeze(SigmaX(5,5,:))), 'g-',...
       t, sqrt(squeeze(SigmaX(6,6,:))), 'b-');
  xlabel('Time');
  ylabel('Velocity Error, per Axis');
  title(strcat(title_hdr, ' Estimated Velocity Error and 95% Error Bounds'));

