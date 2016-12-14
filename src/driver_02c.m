%
% Copyright 2016 Kurt Motekew
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%

% This driver does a sanity check for processing different observation types
% for both the U-D and SRIF filters.  U-D obs are still processed individually
% while SRIF obs are processed both individually and as measurement sets.
% A QR version of the SRIF is also performed.

close all;
clear;

SF95_3D = 2.796;
k2_3d = SF95_3D*SF95_3D;

  % INPUTS
ntest = 20;
  % True location of object to be located
rho = [0.25 0.25 0]';
  % Tracker locations
blen = 1;
tkrs = [
         0       0       blen ;
         blen    0       blen ;
         blen    blen    blen ;
         0       blen    blen
      ]';
nmeas = size(tkrs,2);                        % Batch init
nmeas2 = 100;                                % Sequential obs, static

  % Tracker accuracy
srng = .05;                                      % .05 DU range uncertainty
sang = pi*1/180;                                 % 1 deg angle uncertainty
vrng = srng*srng;
W = 1/vrng;
Wsqrt = 1/srng;
y(nmeas) = 0;
y2(3,nmeas2) = 0;
testnum = 1:ntest;
miss_ud = zeros(1,ntest);
miss_srif = zeros(1,ntest);
miss_srifb = zeros(1,ntest);
miss_qr = zeros(1,ntest);
contained_3d_ud = 0;
contained_3d_srif = 0;
contained_3d_srifb = 0;
contained_3d_qr = 0;
ud_time = 0;
srif_time = 0;
srifb_time = 0;
qr_time = 0;

%
% Loop over the number of trials
%

for jj = 1:ntest
    % Create synthetic measurements using error budget
    % First batch is to initialize the filters
  for ii = 1:nmeas
    s = rho - tkrs(:,ii);
    yerr = srng*randn;
    y(ii) = norm(s) + yerr;
  end
    % Additional observations - range and pointing
  for ii = 1:nmeas2
      % Determine which tracker to get the obs from
    itkr = mod(ii,nmeas);
    if itkr == 0
      itkr = nmeas;
    end
    s = rho - tkrs(:,itkr);
    serr = srng*randn;
    smag = norm(s);
    y2(1,ii) = smag + serr;
    shat = s/smag;
    uerr = sang*randn;
    verr = sang*randn;
    y2(2:3,ii) = shat(1:2,1) + [uerr ; verr];
  end

    % Estimate location using initial set of obs - range only
  [phat0, SigmaP0, ~] = box_locate(tkrs, y, W);
    % Sequential estimation based on initial estimate
  phat_ud = phat0;
  [U, D] = mth_udut2(SigmaP0);                   % U-D, SigmaP = UDU'
    % Use Householder method for initial estimate and get info array
    % Range only for init as with Householder
  [phat_srif, SigmaP_srif, R, z, ~] = box_locate_hh(tkrs, y, Wsqrt*eye(nmeas));
  phat_srifb = phat_srif;
  Rb = R;
  z = 0*z;
    % QR init, similar to Householder
  [phat_qr, SigmaP_qr, Rqr, qty, ~] = box_locate_qr(tkrs, y, Wsqrt*eye(nmeas));
  qty = 0*qty;

  %
  % Updates - use range and pointing obs
  %

    % U-D
  vrpnt = [srng*srng sang*sang sang*sang];       % Range and pointing variances
  tic;
  for ii = 1:nmeas2
    itkr = mod(ii,nmeas);
    if itkr == 0
      itkr = nmeas;
    end
    Ap = est_drpnt_dloc(tkrs(:,itkr), phat_ud);
    s = phat_ud - tkrs(:,itkr);
    smag = norm(s);
    shat = s/smag;
    yc = [smag ; shat(1:2,1)];
    r = y2(:,ii) - yc;
      % Loop through each obs in the set
    for kk = 1:3
      [phat_ud, U, D] = est_upd_ud(phat_ud, U, D, Ap(kk,:), r(kk), vrpnt(kk));
    end
  end
  ud_time = ud_time + toc;

    % SRIF
  Wsqrts = [1/srng 1/sang 1/sang];
  tic;
  for ii = 1:nmeas2
    itkr = mod(ii,nmeas);
    if itkr == 0
      itkr = nmeas;
    end
    Ap = est_drpnt_dloc(tkrs(:,itkr), phat_srif);
    s = phat_srif - tkrs(:,itkr);
    smag = norm(s);
    shat = s/smag;
    yc = [smag ; shat(1:2,1)];
    r = y2(:,ii) - yc;
      % Loop through each obs in the set
    for kk = 1:3
      [dp, R, z, ~] = est_upd_hhsrif(R, z, Ap(kk,:), r(kk), Wsqrts(kk));
      z = 0*z;
      phat_srif = phat_srif + dp;
    end
  end
  srif_time = srif_time + toc;

    % SRIF Batch
  Wsqrtb = [1/srng 0 0 ; 0 1/sang 0 ; 0 0 1/sang];
  tic;
  for ii = 1:nmeas2
    itkr = mod(ii,nmeas);
    if itkr == 0
      itkr = nmeas;
    end
    Ap = est_drpnt_dloc(tkrs(:,itkr), phat_srifb);
    s = phat_srifb - tkrs(:,itkr);
    smag = norm(s);
    shat = s/smag;
    yc = [smag ; shat(1:2,1)];
    r = y2(:,ii) - yc;
      % Process measurement sets
    [dp, Rb, z, ~] = est_upd_hhsrif(Rb, z, Ap, r, Wsqrtb);
    z = 0*z;
    phat_srifb = phat_srifb + dp;
  end
  srifb_time = srifb_time + toc;

    % QR Batch
  tic;
  for ii = 1:nmeas2
    itkr = mod(ii,nmeas);
    if itkr == 0
      itkr = nmeas;
    end
    Ap = est_drpnt_dloc(tkrs(:,itkr), phat_qr);
    s = phat_qr - tkrs(:,itkr);
    smag = norm(s);
    shat = s/smag;
    yc = [smag ; shat(1:2,1)];
    r = y2(:,ii) - yc;
      % Process measurement sets
    [dp, Rqr, qty] = est_upd_qrsrif(Rqr, qty, Ap, r, Wsqrtb);
    qty = 0*qty;
    phat_qr = phat_qr + dp;
  end
  qr_time = qr_time + toc;

  %
  % End updates
  %

    % Get containment stats for each
  miss_ud(jj) = norm(phat_ud - rho);
  if (SF95_3D > mth_mahalanobis(rho, phat_ud, U*D*U'))
    contained_3d_ud = contained_3d_ud + 1;
  end
    %
  miss_srif(jj) = norm(phat_srif - rho);
  Rinv = mth_triinv(R);
  if (SF95_3D > mth_mahalanobis(rho, phat_srif, Rinv*Rinv'))
    contained_3d_srif = contained_3d_srif + 1;
  end
    %
  miss_srifb(jj) = norm(phat_srifb - rho);
  Rbinv = mth_triinv(Rb);
  if (SF95_3D > mth_mahalanobis(rho, phat_srifb, Rbinv*Rbinv'))
    contained_3d_srifb = contained_3d_srifb + 1;
  end
    %
  miss_qr(jj) = norm(phat_qr - rho);
  Rqrinv = mth_triinv(Rqr);
  if (SF95_3D > mth_mahalanobis(rho, phat_qr, Rqrinv*Rqrinv'))
    contained_3d_qr = contained_3d_qr + 1;
  end
end
p95_3d_ud = 100*contained_3d_ud/ntest;
p95_3d_srif = 100*contained_3d_srif/ntest;
p95_3d_srifb = 100*contained_3d_srifb/ntest;
p95_3d_qr = 100*contained_3d_qr/ntest;

figure; hold on;
plot(testnum, miss_ud, 'o', testnum, miss_srif, '+',...
     testnum, miss_srifb, '*', testnum, miss_qr, 'o');
xlabel('Trial');
ylabel('RSS Miss Distance');
legend('U-D', 'SRIF', 'SRIF B', 'QR');

fprintf('\nU-D containment: %1.1f', p95_3d_ud);
fprintf(' in %1.4f seconds', ud_time);
fprintf('\nSRIF containment: %1.1f', p95_3d_srif);
fprintf(' in %1.4f seconds', srif_time);
fprintf('\nSRIF Batch containment: %1.1f', p95_3d_srifb);
fprintf(' in %1.4f seconds', srifb_time);
fprintf('\nQR Batch containment: %1.1f', p95_3d_qr);
fprintf(' in %1.4f seconds', qr_time);

fprintf('\n');


