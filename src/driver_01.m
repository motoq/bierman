%
% Copyright 2016 Kurt Motekew
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%

% This driver script compares full batch, sequential batch, stabilized Kalman,
% Potter mechanization, and U-D sequential estimation methods.  It is calling
% filtering functions that are put together to illustrate their form and make
% comparisons easier.  The driver_02 script calls functions that are more
% gearded towards speed.


close all;
clear;

SF95_3D = 2.796;
k2_3d = SF95_3D*SF95_3D;

  % INPUTS
ntest = 20;
  % Target location
rho = [0.25 0.25 0]';
  % Tracker locations
blen = 1;
tkr = [
        0       0       blen ;
        blen    0       blen ;
        blen    blen    blen ;
        0       blen    blen
      ]';
nmeas = size(tkr,2);                        % Batch init
nmeas2 = 100;                                % Sequential obs, static
tkr2 = zeros(3,nmeas2);

% Plot geometry
box_plot(rho, tkr, blen);
view([70 20]);

  % Tracker accuracy
srng = .05;
vrng = srng*srng;
W = vrng^-1;         
y(nmeas) = 0;
y2(nmeas2) = 0;
YERRS(nmeas+nmeas2,ntest) = 0;
locs0 = zeros(3,ntest);
testnum = 1:ntest;
miss0 = zeros(1,ntest);
miss_sb = zeros(1,ntest);
miss_fb = zeros(1,ntest);
miss_kf = zeros(1,ntest);
miss_pt = zeros(1,ntest);
miss_ud = zeros(1,ntest);
contained_3d0 = 0;
contained_3d_sb = 0;
contained_3d_fb = 0;
contained_3d_kf = 0;
contained_3d_pt = 0;
contained_3d_ud = 0;

%
% Loop over the number of trials
%

for jj = 1:ntest
    % Create synthetic measurements using error budget
  for ii = 1:nmeas
    s = rho - tkr(:,ii);
    yerr = srng*randn;
    y(ii) = norm(s) + yerr;
    YERRS(ii,jj) = abs(yerr);
  end
    % Additional observations
  for ii = 1:nmeas2
      % Determine which tracker to get the obs from
    itkr = mod(ii,nmeas);
    if itkr == 0
      itkr = nmeas;
    end
    s = rho - tkr(:,itkr);
    yerr = srng*randn;
    y2(ii) = norm(s) + yerr;
    YERRS(ii+nmeas,jj) = abs(yerr);
  end
    % Estimate location using initial set of obs
  [phat0, SigmaP0, ~] = box_locate(tkr, y, W);
  locs0(:,jj) = phat0;
  if (SF95_3D > mth_mahalanobis(rho, phat0, SigmaP0))
    contained_3d0 = contained_3d0 + 1;
  end
  miss0(jj) = norm(phat0 - rho);

    % Sequential estimation based on initial estimate
  SigmaP_sb = SigmaP0;                           % Sequential batch
  phat_sb = phat0;
  SigmaP_kf = SigmaP0;                           % Kalman filter
  phat_kf = phat0;
  SigmaP_pt = SigmaP0;                           % Potter form
  phat_pt = phat0;
  SigmaP_ud = SigmaP0;                           % Potter form
  phat_ud = phat0;
  for ii = 1:nmeas2
      % Determine which tracker to get the obs from
    itkr = mod(ii,nmeas);
    if itkr == 0
      itkr = nmeas;
    end
    tkr2(:,ii) = tkr(:,itkr);
    [phat_sb, SigmaP_sb] = box_update(phat_sb, SigmaP_sb,...
                                      tkr2(:,ii), y2(ii), W,...
                                      'SEQBATCH');
    [phat_kf, SigmaP_kf] = box_update(phat_kf, SigmaP_kf,...
                                      tkr2(:,ii), y2(ii), W,...
                                      'KALMAN');
    [phat_pt, SigmaP_pt] = box_update(phat_pt, SigmaP_pt,...
                                      tkr2(:,ii), y2(ii), W,...
                                      'POTTER');
    [phat_ud, SigmaP_ud] = box_update(phat_ud, SigmaP_ud,...
                                      tkr2(:,ii), y2(ii), W, 'UD');
  end

    % Get containment stats for each
  if (SF95_3D > mth_mahalanobis(rho, phat_sb, SigmaP_sb))
    contained_3d_sb = contained_3d_sb + 1;
  end
  miss_sb(jj) = norm(phat_sb - rho);
  if (SF95_3D > mth_mahalanobis(rho, phat_kf, SigmaP_kf))
    contained_3d_kf = contained_3d_kf + 1;
  end
  miss_kf(jj) = norm(phat_kf - rho);
  if (SF95_3D > mth_mahalanobis(rho, phat_pt, SigmaP_pt))
    contained_3d_pt = contained_3d_pt + 1;
  end
  miss_pt(jj) = norm(phat_pt - rho);
  if (SF95_3D > mth_mahalanobis(rho, phat_ud, SigmaP_ud))
    contained_3d_ud = contained_3d_ud + 1;
  end
  miss_ud(jj) = norm(phat_ud - rho);

  [phat_fb, SigmaP_fb, ~] = box_locate([tkr tkr2], [y y2], W);
  if (SF95_3D > mth_mahalanobis(rho, phat_fb, SigmaP_fb))
    contained_3d_fb = contained_3d_fb + 1;
  end
  miss_fb(jj) = norm(phat_fb - rho);
end
p95_3d0 = 100*contained_3d0/ntest;
p95_3d_sb = 100*contained_3d_sb/ntest;
p95_3d_fb = 100*contained_3d_fb/ntest;
p95_3d_kf = 100*contained_3d_kf/ntest;
p95_3d_pt = 100*contained_3d_pt/ntest;
p95_3d_ud = 100*contained_3d_ud/ntest;
p95_3d0s = sprintf('%4.1f%%', 100*contained_3d0/ntest);

  % 3D covariance
ApTWAp = box_infom(tkr, rho, W);
SigmaP = ApTWAp^-1;
  % 95% Covariance
SigmaP3D = k2_3d*SigmaP;

figure; hold on;
for jj=1:ntest
  scatter3(locs0(1,jj), locs0(2,jj), locs0(3,jj), 'r');
end
[XX, YY, ZZ] = matrix3X3_points(SigmaP3D, 20);
mesh(XX + rho(1,1), YY + rho(2,1), ZZ + rho(3,1));
hidden('off');
xlabel('X');
ylabel('Y');
zlabel('Z');
tstring = strcat('3D Estimates, 95% Covariance:  ', p95_3d0s, ' Containment');
title(tstring);
view([70 20]);

figure; hold on;
plot(testnum, miss0, '-x', testnum, miss_sb, '-*', testnum, miss_fb, 'o',...
     testnum, miss_kf, 'd', testnum, miss_pt, 'd', testnum, miss_ud, 's');
xlabel('Trial');
ylabel('RSS Miss Distance');
legend('Init', 'Seq Batch', 'Full Batch', 'Kalman', 'Potter', 'U-D');

figure; hold on;
imagesc(YERRS);
xlabel('Trial');
ylabel('Obs#');
title('Range Error');
colorbar;
%axis image;

fprintf('\nBatch initializer containment: %1.1f', p95_3d0);
fprintf('\nSequential batch containment: %1.1f', p95_3d_sb);
fprintf('\nKalman containment: %1.1f', p95_3d_kf);
fprintf('\nPotter containment: %1.1f', p95_3d_pt);
fprintf('\nU-D containment: %1.1f', p95_3d_ud);
fprintf('\nFull batch containment: %1.1f', p95_3d_fb);

fprintf('\n');


