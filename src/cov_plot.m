function cov_plot(title_hdr, t, SigmaX)
% COV_PLOT plots position and velocity 95% covariance.
%
%-----------------------------------------------------------------------
% Copyright 2018 Kurt Motekew
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%-----------------------------------------------------------------------
%
% Inputs
%   title_hdr   String indicating what method is used
%   t           Array of times
%   SigmaX      State covariance, [6x6xN]
%
% Author:  Kurt Motekew    20180321
%

    % Scale covariance by 95% containment
  SF95_3D = 2.796;
  SigmaX = SF95_3D*SF95_3D*SigmaX;

    % Position
  figure; hold on;
  plot(t, sqrt(squeeze(SigmaX(1,1,:))), 'r-',...
       t, sqrt(squeeze(SigmaX(2,2,:))), 'g-',...
       t, sqrt(squeeze(SigmaX(3,3,:))), 'b-');
  legend('x', 'y', 'z');
  xlabel('Time');
  ylabel('Position Error, per Axis');
  title(strcat(title_hdr, ' 95% Position Error Bounds'));

    % Velocity
  figure; hold on;
  plot(t, sqrt(squeeze(SigmaX(4,4,:))), 'r-',...
       t, sqrt(squeeze(SigmaX(5,5,:))), 'g-',...
       t, sqrt(squeeze(SigmaX(6,6,:))), 'b-');
  legend('x', 'y', 'z');
  xlabel('Time');
  ylabel('Velocity Error, per Axis');
  title(strcat(title_hdr, ' 95% Velocity Error Bounds'));

