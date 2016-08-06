function SigmaP3D = box_plot(rho, r, blen)
% BOX_PLOT Plots tracker location and an exagerated covariance around
% the tracked object location for visualization purposes.  This function
% creates a new figure.  Everything is expected to be in the first quadrant.
%
%-----------------------------------------------------------------------
% Copyright 2016 Kurt Motekew
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%-----------------------------------------------------------------------
%
% Inputs:
%  rho    Tracked object location, [3x1]
%  r      [3xN] list of N tracker locations
%  blen   Maximum length of the boxed area the trackers reside
%
% Kurt Motekew   2016/06/22
%

nmeas = size(r,2);
srng = blen/10;

  % Tracker accuracy for plotting
SigmaRng = srng*srng;
W = SigmaRng^-1;

  % 3D solution
ApTWAp = zeros(3);
for ii = 1:nmeas
  Api = drng_dloc(r(:,ii), rho);
  ApTWAp = ApTWAp + Api'*W*Api;
end
SigmaP3D = ApTWAp^-1;

figure; hold on;
axis([0 blen 0 blen]);
[XX, YY, ZZ] = matrix3X3_points(SigmaP3D, 20);
mesh(XX + rho(1,1), YY + rho(2,1), ZZ + rho(3,1));
hidden('off');

r2 = srng/2;
r2 = r2*r2;
TrackerModel = [r2  0  0 ;
                 0 r2  0 ;
                 0  0 r2 ];
[XX, YY, ZZ] = matrix3X3_points(TrackerModel, 20);
for ii = 1:nmeas
  surf(XX + r(1,ii), YY + r(2,ii), ZZ + r(3,ii));
end
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Boxed Tracker Geometry');
box on;
