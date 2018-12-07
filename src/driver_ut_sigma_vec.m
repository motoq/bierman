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

alpha = .69;                           % Small positive value, 1e-4 <= a <= 1
kappa = 0; %3 - L;                     % Secondary scaling parameter
beta  = 2;

[Chi, w_m, w_c] = est_ut_sigma_vec(phat0, SigmaP0, alpha, kappa, beta);
n_params = size(Chi, 1);
n_sigma_vec = size(Chi, 2);

  % Plot covariance and sigma vectors
figure; hold on;
matrix3X3_plot(SigmaP0, 20, false);
ChiPlt = Chi';
scatter3(ChiPlt(:,1), ChiPlt(:,2), ChiPlt(:,3), 'r', 'filled');
title('Unscented Transformation Sigma Vectors');

