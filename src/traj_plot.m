function traj_plot(rhom, rho, r, blen)
% TRAJ_PLOT plots tracker locations and the trajectory of a moving
% tracked object ("ideal" and "truth" trajectories) to illustrate the
% the scenario being modeled.  The ideal trajectory is yellow, the true
% trajectory is blue, and the common starting location is a red point.
%
%-----------------------------------------------------------------------
% Copyright 2014 Kurt Motekew
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%-----------------------------------------------------------------------
%
% Inputs:
%  rhom   [3xN1] list of the tracked object ideal trajectory points
%  rho    [3xN2] list of the tracked object true trajectory points
%  r      [3xN3] list of N3 tracker locations
%  blen   Maximum length (size) of the boxed area the trackers reside
%
% Kurt Motekew   2014/11/05
%

  % Plot tracked object locations - mark first location as red
figure; hold on;
  % Ideal
ntimes = size(rhom,2);
for ii = 1:ntimes
  scatter3(rhom(1,ii), rhom(2,ii), rhom(3,ii), 'y');
end
  % Actual
ntimes = size(rho,2);
for ii = 1:ntimes
  scatter3(rho(1,ii), rho(2,ii), rho(3,ii), 'b');
end
scatter3(rho(1,1), rho(2,1), rho(3,1), 'r')

  % Plot tracker locations
ntkr = size(r,2);
srng = blen/10;

r2 = srng/2; 
r2 = r2*r2;
TrackerModel = [r2  0  0 ;
                 0 r2  0 ;
                 0  0 r2 ];
[XX, YY, ZZ] = matrix3X3_points(TrackerModel, 20);
for ii = 1:ntkr
  surf(XX + r(1,ii), YY + r(2,ii), ZZ + r(3,ii));
end
axis([0 blen 0 blen]);
xlabel('X');
ylabel('Y');
zlabel('Z');
