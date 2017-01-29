%
% Copyright 2016 Kurt Motekew
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%

% This driver script does speed comparisions with the Kalman, U-D, and SRIF
% estimation methods.  Containment is also output


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
SigmaRng = srng*srng;
Wsqrt = 1/srng;
y(nmeas) = 0;
y2(nmeas2) = 0;
testnum = 1:ntest;
miss_srif = zeros(1,ntest);
miss_srifb = zeros(1,ntest);
miss_srif_sch = zeros(1,ntest);
miss_srif_batch = zeros(1,ntest);
contained_3d_srif = 0;
contained_3d_srifb = 0;
contained_3d_srif_sch = 0;
contained_3d_srif_batch = 0;
srif_time = 0;
srifb_time = 0;
srif_sch_time = 0;
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

    % SRIF
  phat_srif = phat_srif0;
  R = R0;
  z = z0;
  tic;
  for ii = 1:nmeas2
    itkr = mod(ii,nmeas);
    if itkr == 0
      itkr = nmeas;
    end
    Ap = est_drng_dloc(tkrs(:,itkr), phat_srif);
    yc = norm(phat_srif - tkrs(:,itkr));
    r = y2(ii) - yc;
    [dp, R, z, ~] = est_upd_hhsrif(R, z, Ap, r, Wsqrt);
    z = 0*z;
    phat_srif = phat_srif + dp;
  end
  srif_time = srif_time + toc;

    % SRIF with Bias
  phat_srifb = phat_srif0;
  Rx = R0;
  Rxy = zeros(3,3);                 % 3 solve for, 3 pos bias per obs
  z = z0;
  Ry = sqrtm((SigmaBlen*eye(3))^-1);
  zy = zeros(3,1);
  tic;
  for ii = 1:nmeas2
    itkr = mod(ii,nmeas);
    if itkr == 0
      itkr = nmeas;
    end
    Ap = est_drng_dloc(tkrs(:,itkr), phat_srifb);
    yc = norm(phat_srifb - tkrs(:,itkr));
    r = y2(ii) - yc;
    Ay = est_drng_dpos(tkrs(:,itkr), phat_srifb);
    [dp, Rx, Rxy, z, Ry, zy] = est_upd_hhsrif_bias(Rx, Rxy, z, Ap, r, Wsqrt,...
                                                               Ry, zy, Ay);
    z = 0*z;
    phat_srifb = phat_srifb + dp;
  end
  srifb_time = srifb_time + toc;

    % SRIF with Schmidt mapping
  phat_srif_sch = phat_srif0;
  Rsch = R0;
  z = z0;
  tic;
  for ii = 1:nmeas2
    itkr = mod(ii,nmeas);
    if itkr == 0
      itkr = nmeas;
    end
    Ap = est_drng_dloc(tkrs(:,itkr), phat_srif_sch);
    yc = norm(phat_srif_sch - tkrs(:,itkr));
    r = y2(ii) - yc;
    Aq = est_drng_dpos(tkrs(:,itkr), phat_srif_sch);
    Wsqrt = sqrtm((SigmaRng + Aq*SigmaBlen*Aq')^-1);
    [dp, Rsch, z, ~] = est_upd_hhsrif(Rsch, z, Ap, r, Wsqrt);
    z = 0*z;
    phat_srif_sch = phat_srif_sch + dp;
  end
  srif_sch_time = srif_sch_time + toc;

    % SRIF full batch
  tic;
  [phat_batch, SigmaP_batch, ~, ~, ~] = box_locate_hh(tkr_pos_all, y_all,...
                                            Wsqrt*eye(nmeas+nmeas2), sblen);
  srif_batch_time = srif_batch_time + toc;


  %
  % End updates
  %

    % Get containment stats for each
  miss_srif(jj) = norm(phat_srif - rho);
  Rinv = mth_triinv(R);
  if (SF95_3D > mth_mahalanobis(rho, phat_srif, Rinv*Rinv'))
    contained_3d_srif = contained_3d_srif + 1;
  end
  miss_srifb(jj) = norm(phat_srifb - rho);
  Rn = [ Rx Rxy ; zeros(3) Ry ];
  Rninv = mth_triinv(Rn);
  Pn = Rninv*Rninv';
  SigmaX = Pn(1:3,1:3);
  %Pxy =  Pn(1:3,4:6);
  Pny = Pn(4:6,4:6);
  SigmaX = SigmaX + Pny;  % How to map for non-identity obs/solve for
  if (SF95_3D > mth_mahalanobis(rho, phat_srifb, SigmaX))
    contained_3d_srifb = contained_3d_srifb + 1;
  end
  miss_srif_sch(jj) = norm(phat_srif_sch - rho);
  Rinv = mth_triinv(Rsch);
  if (SF95_3D > mth_mahalanobis(rho, phat_srif_sch, Rinv*Rinv'))
    contained_3d_srif_sch = contained_3d_srif_sch + 1;
  end
  miss_srif_batch(jj) = norm(phat_batch - rho);
  if (SF95_3D > mth_mahalanobis(rho, phat_batch, SigmaP_batch))
    contained_3d_srif_batch = contained_3d_srif_batch + 1;
  end
end
p95_3d_srif = 100*contained_3d_srif/ntest;
p95_3d_srifb = 100*contained_3d_srifb/ntest;
p95_3d_srif_sch = 100*contained_3d_srif_sch/ntest;
p95_3d_srif_batch = 100*contained_3d_srif_batch/ntest;

figure; hold on;
plot(testnum, miss_srif, 's', testnum, miss_srifb, 'x',...
     testnum, miss_srif_sch, '*', testnum, miss_srif_batch, 'o');
xlabel('Trial');
ylabel('RSS Miss Distance');
legend('SRIF', 'SRIF B', 'Schmidt SRIF', 'Full Batch');

fprintf('\nSRIF containment: %1.1f', p95_3d_srif);
fprintf(' in %1.4f seconds', srif_time);
fprintf('\nSRIF B containment: %1.1f', p95_3d_srifb);
fprintf(' in %1.4f seconds', srifb_time);
fprintf('\nSch SRIF containment: %1.1f', p95_3d_srif_sch);
fprintf(' in %1.4f seconds', srif_sch_time);
fprintf('\nFull Batch SRIF containment: %1.1f', p95_3d_srif_batch);
fprintf(' in %1.4f seconds', srif_batch_time);

fprintf('\n');


