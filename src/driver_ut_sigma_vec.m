%
% Copyright 2018 Kurt Motekew
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%

%
% This driver script illustrates the unscented transformation of a 3D estimate
% given an associated covariance.
%

clear;
close all;

  % Estimate (origin of reference frame) and covariance
phat0 = [ 0 0 0]';
SigmaP0 = [ 3 2 1 ; 2 4 3 ; 1 3 5 ];

L = size(phat0,1);                     % Dimension of state vector
Chi = zeros(L, 2*L +1);                % Sigma vectors
Chi(:,1) = phat0;

alpha = .69;                           % Small positive value, 1 <= a <= 1e-4
kappa = 0; %3 - L;                     % Secondary scaling parameter
lambda = alpha*alpha*(L + kappa) - L;

  % Scale estimate covariance and form sigma vectors
SigmaP0scaled = sqrtm((L + lambda)*SigmaP0);
for ii = 2:(L+1)
  Chi(:,ii) = phat0 + SigmaP0scaled(:,ii-1);
  Chi(:,ii+L) = phat0 - SigmaP0scaled(:,ii-1);
end

  % Plot covariance and sigma vectors
figure; hold on;
matrix3X3_plot(SigmaP0, 20, false);
ChiPlt = Chi';
scatter3(ChiPlt(:,1), ChiPlt(:,2), ChiPlt(:,3), 'r', 'filled');
title('Unscented Transformation Sigma Vectors');
