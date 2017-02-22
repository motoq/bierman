%
% Copyright 2017 Kurt Motekew
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%

% This driver adds bias to the measurement model by perturbing the sensor
% locations.  It then performs accuracy, containment, and speed comparisons
% of the Schmodt Consider Kalman(SCK), SRIF, SRIF with 'a posteriori'
% covariance inflation, SRIF with bias inclusion during processing, and SRIF
% full batch processing.  Speed and containment are compared.  Key comparisions
% are the lack of containment with no bias (SRIF), speed and containment 
% between SCK and SRIF B, and full batch SRIF vs. SRIF 'a posteriori' 
% containment.  While the full batch method does not improve the estimate
% for this scenario (it can often provide an order of magnitude improvement
% in the miss distance over other methods) it does have the most reliable
% estimate covariance.


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
sblen = .003*blen;                          % .3% of total distance for tracker
SigmaBlen = sblen*sblen;
SigmaY = SigmaBlen*eye(3);
tkrs = [                                    % uncertainty
         0       0       blen ;
         blen    0       blen ;
         blen    blen    blen ;
         0       blen    blen
      ]';
nmeas = size(tkrs,2);                        % Batch init
nmeas2 = 20;                                % Sequential obs, static

  % Tracker accuracy
srng = .001;
SigmaZ = srng*srng;
Wsqrt = 1/srng;
W = Wsqrt*Wsqrt';
y(nmeas) = 0;
y2(nmeas2) = 0;
testnum = 1:ntest;
miss_schmidt = zeros(1,ntest);
miss_srif = zeros(1,ntest);
miss_srif_apost = zeros(1,ntest);
miss_srifb = zeros(1,ntest);
miss_srif_batch = zeros(1,ntest);
contained_3d_schmidt = 0;
contained_3d_srif = 0;
contained_3d_srif_apost = 0;
contained_3d_srifb = 0;
contained_3d_srif_batch = 0;
schmidt_time = 0;
srif_time = 0;
srif_apost_time = 0;
srifb_time = 0;
srif_batch_time = 0;

tkr_pos_all = zeros(3,nmeas+nmeas2);

%
% Loop over the number of trials
%

for jj = 1:ntest
  tkr_deltas = sblen*randn(3,nmeas);
    % Create synthetic measurements using error budget
  for ii = 1:nmeas
    tkr_pos_all(:,ii) = tkrs(:,ii);
    s = rho - tkrs(:,ii) + tkr_deltas(:,ii);
    yerr = srng*randn;
    y(ii) = norm(s) + yerr;
  end
    % Additional observations
  for ii = 1:nmeas2
      % Determine which tracker to get the obs from
    itkr = mod(ii,nmeas);
    if itkr == 0
      itkr = nmeas;
    end
    tkr_pos_all(:,nmeas+ii) = tkrs(:,itkr);
    s = rho - tkrs(:,itkr) + tkr_deltas(:,itkr);
    yerr = srng*randn;
    y2(ii) = norm(s) + yerr;
  end
  y_all(nmeas+1:nmeas+nmeas2) = y2;
  y_all(1:nmeas) = y;
    % Use Householder method for initial estimate and get info array
  [phat_srif0, SigmaP_srif0, R0, z0, ~] =...
                                      box_locate_hh(tkrs, y,...
                                                    Wsqrt*eye(nmeas), sblen);
  z0 = 0*z0;

  %
  % Updates
  %

    % Schmidt consider Kalman
  phat_schmidt = phat_srif0;
  SigmaP_schmidt = SigmaP_srif0;
  SigmaPY = zeros(3,3);
  tic;
  for ii = 1:nmeas2
    itkr = mod(ii,nmeas);
    if itkr == 0
      itkr = nmeas;
    end
    Ax = est_drng_dloc(tkrs(:,itkr), phat_schmidt);
    Ay = est_drng_dpos(tkrs(:,itkr), phat_schmidt);
    yc = norm(phat_schmidt - tkrs(:,itkr));
    r = y2(ii) - yc;
    [phat_schmidt, SigmaP_schmidt, SigmaPY] = est_upd_schmidt(...
                                        phat_schmidt, SigmaP_schmidt, Ax,...
                                        SigmaY, Ay, SigmaPY, r, SigmaZ);
  end
  schmidt_time = schmidt_time + toc;

    % SRIF with no bias incorporation
  phat_srif = phat_srif0;
  R = R0;
  z = z0;
  tic;
  for ii = 1:nmeas2
    itkr = mod(ii,nmeas);
    if itkr == 0
      itkr = nmeas;
    end
    Ax = est_drng_dloc(tkrs(:,itkr), phat_srif);
    yc = norm(phat_srif - tkrs(:,itkr));
    r = y2(ii) - yc;
    [dp, R, z, ~] = est_upd_hhsrif(R, z, Ax, r, Wsqrt);
    z = 0*z;
    phat_srif = phat_srif + dp;
  end
  srif_time = srif_time + toc;

    % SRIF with 'a posteriori' bias inflation
  phat_srif_apost = phat_srif0;
  Rapost = R0;
  z = z0;
  AxTWAy = zeros(3); % apost
  tic;
  for ii = 1:nmeas2
    itkr = mod(ii,nmeas);
    if itkr == 0
      itkr = nmeas;
    end
    Ax = est_drng_dloc(tkrs(:,itkr), phat_srif_apost);
    Ay = est_drng_dpos(tkrs(:,itkr), phat_srif_apost);
    AxTWAy = AxTWAy + Ax'*W*Ay;
    yc = norm(phat_srif_apost - tkrs(:,itkr));
    r = y2(ii) - yc;
    [dp, Rapost, z, ~] = est_upd_hhsrif(Rapost, z, Ax, r, Wsqrt);
    z = 0*z;
    phat_srif_apost = phat_srif_apost + dp;
  end
  srif_apost_time = srif_apost_time + toc;

    % SRIF with Bias
  phat_srifb = phat_srif0;
  Rx = R0;
  Rxy = zeros(3,3);                 % 3 solve for, 3 pos bias per obs
  z = z0;
  Ry = sqrtm(SigmaY^-1);
  zy = zeros(3,1);
  tic;
  for ii = 1:nmeas2
    itkr = mod(ii,nmeas);
    if itkr == 0
      itkr = nmeas;
    end
    Ax = est_drng_dloc(tkrs(:,itkr), phat_srifb);
    Ay = est_drng_dpos(tkrs(:,itkr), phat_srifb);
    yc = norm(phat_srifb - tkrs(:,itkr));
    r = y2(ii) - yc;
    [dp, Rx, Rxy, z, Ry, zy] = est_upd_hhsrif_bias(Rx, Rxy, z, Ax, r, Wsqrt,...
                                                               Ry, zy, Ay);
    z = 0*z;
    zy = 0*zy;
    phat_srifb = phat_srifb + dp;
  end
  srifb_time = srifb_time + toc;

    % SRIF full batch
  tic;
  [phat_batch, SigmaP_batch, ~, ~, ~] = box_locate_hh(tkr_pos_all, y_all,...
                                            Wsqrt*eye(nmeas+nmeas2), sblen);
  srif_batch_time = srif_batch_time + toc;


  %
  % End updates
  %

    % Get containment stats for each
  miss_schmidt(jj) = norm(phat_schmidt - rho);
  if (SF95_3D > mth_mahalanobis(rho, phat_schmidt, SigmaP_schmidt))
    contained_3d_schmidt = contained_3d_schmidt + 1;
  end
    %
  miss_srif(jj) = norm(phat_srif - rho);
  Rinv = mth_triinv(R);
  SigmaX = Rinv*Rinv';
  if (SF95_3D > mth_mahalanobis(rho, phat_srif, SigmaX))
    contained_3d_srif = contained_3d_srif + 1;
  end
    %
  miss_srif_apost(jj) = norm(phat_srif_apost - rho);
  Rinv = mth_triinv(Rapost);
  SigmaX = Rinv*Rinv';
  SigmaX = SigmaX + SigmaX*AxTWAy*SigmaY*AxTWAy'*SigmaX;
  if (SF95_3D > mth_mahalanobis(rho, phat_srif_apost, SigmaX))
    contained_3d_srif_apost = contained_3d_srif_apost + 1;
  end
    %
  miss_srifb(jj) = norm(phat_srifb - rho);
  Rxinv = mth_triinv(Rx);
  S = -Rxinv*Rxy;
  Phatc = Rxinv*Rxinv';
  Ryinv = mth_triinv(Ry);
  Py = Ryinv*Ryinv';
  SigmaX = Phatc + S*Py*S';
  if (SF95_3D > mth_mahalanobis(rho, phat_srifb, SigmaX))
    contained_3d_srifb = contained_3d_srifb + 1;
  end
  miss_srif_batch(jj) = norm(phat_batch - rho);
  if (SF95_3D > mth_mahalanobis(rho, phat_batch, SigmaP_batch))
    contained_3d_srif_batch = contained_3d_srif_batch + 1;
  end
end
p95_3d_schmidt = 100*contained_3d_schmidt/ntest;
p95_3d_srif = 100*contained_3d_srif/ntest;
p95_3d_srif_apost = 100*contained_3d_srif_apost/ntest;
p95_3d_srifb = 100*contained_3d_srifb/ntest;
p95_3d_srif_batch = 100*contained_3d_srif_batch/ntest;

figure; hold on;
plot(testnum, miss_schmidt, '+',...
     testnum, miss_srif, 's', testnum, miss_srif_apost, '*',...
     testnum, miss_srifb, 'x', testnum, miss_srif_batch, 'o');
xlabel('Trial');
ylabel('RSS Miss Distance');
legend('SKF', 'SRIF', 'SRIF Apost', 'SRIF B', 'Full Batch');

fprintf('\nSKF containment: %1.1f', p95_3d_schmidt);
fprintf(' in %1.4f seconds', schmidt_time);
fprintf('\nSRIF containment: %1.1f', p95_3d_srif);
fprintf(' in %1.4f seconds', srif_time);
fprintf('\nSRIF a posteriori containment: %1.1f', p95_3d_srif_apost);
fprintf(' in %1.4f seconds', srif_apost_time);
fprintf('\nSRIF B containment: %1.1f', p95_3d_srifb);
fprintf(' in %1.4f seconds', srifb_time);
fprintf('\nFull Batch SRIF containment: %1.1f', p95_3d_srif_batch);
fprintf(' in %1.4f seconds', srif_batch_time);

fprintf('\n');


