%
% Copyright 2016 Kurt Motekew
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%

% This driver script compares compares A WLS solution of the normal equations
% with one using Householder triangularization.  It is primarily meant to
% be a sanity check of the Householder decomposition, using a batch solution
% as an example (vs. a filtered example).


close all;
clear;

SF95_3D = 2.796;
k2_3d = SF95_3D*SF95_3D;

  % INPUTS
ntest = 200;
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
ntkrs = size(tkrs,2);
nmeas = 10;
tkrs2 = zeros(3,nmeas);

% Plot geometry
box_plot(rho, tkrs, blen);
view([70 20]);

  % Tracker accuracy
srng = .05;
vrng = srng*srng;
SigmaO = vrng*eye(nmeas);
W = eye(nmeas)/vrng;
SqrtW = sqrtm(W);

y(nmeas) = 0;
testnum = 1:ntest;
miss_fw = zeros(1,ntest);
miss_qr = zeros(1,ntest);
miss_hh = zeros(1,ntest);
contained_3d_fw = 0;
contained_3d_qr = 0;
contained_3d_hh = 0;
fw_time = 0;
qr_time = 0;
hh_time = 0;

%
% Loop over the number of trials
%

for jj = 1:ntest
    % Create synthetic measurements using error budget
  for ii = 1:nmeas
      % Determine which tracker to get the obs from
    itkr = mod(ii,ntkrs);
    if itkr == 0
      itkr = ntkrs;
    end
    tkrs2(:,ii) = tkrs(:,itkr);
    s = rho - tkrs2(:,ii);
    yerr = srng*randn;
    y(ii) = norm(s) + yerr;
  end
    % Normal equations solution
  tic;
    [phat_fw, SigmaP_fw, ~] = box_locate_fw(tkrs2, y, W);
  fw_time = fw_time + toc;
  if (SF95_3D > mth_mahalanobis(rho, phat_fw, SigmaP_fw))
    contained_3d_fw = contained_3d_fw + 1;
  end
  miss_fw(jj) = norm(phat_fw - rho);
    %
  tic
    [phat_qr, SigmaP_qr, ~, ~, ~] = box_locate_qr(tkrs2, y, SqrtW);
  qr_time = qr_time + toc;
  if (SF95_3D > mth_mahalanobis(rho, phat_qr, SigmaP_qr))
    contained_3d_qr = contained_3d_qr + 1;
  end
  miss_qr(jj) = norm(phat_qr - rho);
    %
  tic
    [phat_hh, SigmaP_hh, ~, ~, ~] = box_locate_hh(tkrs2, y, SqrtW);
  hh_time = hh_time + toc;
  if (SF95_3D > mth_mahalanobis(rho, phat_hh, SigmaP_hh))
    contained_3d_hh = contained_3d_hh + 1;
  end
  miss_hh(jj) = norm(phat_hh - rho);
end
p95_3d_fw = 100*contained_3d_fw/ntest;
p95_3d_qr = 100*contained_3d_qr/ntest;
p95_3d_hh = 100*contained_3d_hh/ntest;

figure; hold on;
plot(testnum, miss_fw, 'o', testnum, miss_qr, 's', testnum, miss_hh, '*');
xlabel('Trial');
ylabel('RSS Miss Distance');
legend('Full Batch', 'QR', 'Householder');

fprintf('\nFull batch containment: %1.1f', p95_3d_fw);
fprintf(' in %1.4f seconds', fw_time);
fprintf('\nQr containment: %1.1f', p95_3d_hh);
fprintf(' in %1.4f seconds', qr_time);
fprintf('\nHouseholder containment: %1.1f', p95_3d_hh);
fprintf(' in %1.4f seconds', hh_time);

fprintf('\n');


