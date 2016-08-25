function [phat, SigmaP] = box_update(phat0, P0, tkr_pos, y, W, mode)
% BOX_UPDATE Updates the geolocation estimate of a stationary object within a
% boxed volume given an a priori estimate, covariance, range tracker location,
% and an additional observation with weighting.  Multiple sequential update
% options are offered.  Each accepts one new scalar observation.  This function
% is written primarily with the intent of illustrating the different update
% methods.  Individual functions exist for selected filtering algorithms that
% are more suited for real world use.
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
%   phat0     A priori position estimate, [3x1]
%   P0        A priori estimate covariance, [3x3]
%   tkr_pos   A [3x1] vector representing a tracker location
%   y         New distance measurements, scalar
%   W         Range uncertianty weighting matrix, 1/sigma_range^2, scalar
%   mode      String indicating what mode to use
%               'SEQBATCH'  Sequential batch mode
%               'KALMAN'    Kalman Filter
%               'POTTER'    Potter Mechanization 
%               'U-D'       U-D
%
% Return:
%   phat     Estimated location, [3x1]
%   SigmaP   Location covariance, [3x3]
%
% Kurt Motekew   2016/07/21
%
    %
    % Generic parameters
    %
    % Moving these three computations externally and passing them instead
    % of tkr_pos and y makes the filtering algorithms below agnostic to the
    % problem as long as iteration is not performed (as long as updating the
    % partials and residual is not necessary).  They are included here for
    % clarity:
    %   1:  Illustrates differential form of linear update
    %   2:  Illustrates (in one place) scaling of observation (residual) and
    %       partials when using a method that assumes the observation
    %       covariance is the identity matrix (normalized form).
    %
    % Ref:  G. J. Bierman, Factorization Methods for
    %       Discrete Sequential Estimation, Dover Publications, Inc.,
    %       Mineola, NY, 1977.
    %
  yc = norm(phat0 - tkr_pos);                    % Computed observation, scalar
  Ap = est_drng_dloc(tkr_pos, phat0);            % Partials, [1x3]
  r = y - yc;                                    % Predicted residual
    % Normalize, Covariance obs = I
  sw = sqrt(W);
  delta = sw*r;                                  % Scaled predicted residual
  Ap = sw*Ap;                                    % Scaled partials
  
    %
    % Sequential update method
    %
  if strcmp('SEQBATCH', mode)
    F0 = P0^-1;                                  % Seqential uses info matrix
    ApTWr =  Ap'*delta;                          % Form normal equations
    ApTWAp = F0 + Ap'*Ap;                        % with information matrix
    SigmaP = ApTWAp^-1;                          % Updated estimate covariance
    dp = SigmaP*ApTWr;                           %dp = SigmaP*(F0*dp + ApTWr);
    phat = phat0 + dp;                           % Updated estimate
  elseif strcmp('KALMAN', mode)
    v = P0*Ap';
    sigma = Ap*v +1;                             % Predicted residual covariance
    K = v/sigma;                                 % Kalman gain
    phat = phat0 + K*delta;                      % State update
    Pbar = P0 - K*v';                            % Optimal covariance update
    v = Pbar*Ap';
    SigmaP = (Pbar - v*K') + K*K';               % Stabilized covariance update
  elseif strcmp('POTTER', mode)
    [S0, ~] = mth_sqrtm(P0);                     % Matrix square root
    vtrans = Ap*S0;                              % [1x3]
    sigma = 1/(vtrans*vtrans' +1);               % Predicted residual covariance
    K = S0*vtrans';                              % Kalman gain
    phat = phat0 + K*(delta*sigma);              % State update
    lambda = sigma/(1 + sqrt(sigma));
    S = S0 - (lambda*K)*vtrans;
    SigmaP = S*S';                               % Stabilized covariance update
  elseif strcmp('UD', mode)
    n = size(phat0,1);
    a = est_drng_dloc(tkr_pos, phat0);           % Not using normalized obs
    sig2 = 1/W;                                  % Variance
    U = mth_udut(P0);                            % Combined U-D
    b = zeros(1,n);
    for jj = n:-1:2
      for kk = 1:(jj-1)
        a(jj) = a(jj) + U(kk,jj)*a(kk);            % a = U'a
      end
      b(jj) = U(jj,jj)*a(jj);                      % b = DU'a
    end
    b(1) = U(1,1)*a(1);

    alpha = sig2 + b(1)*a(1);
    gamma = 1/alpha;
    U(1,1) = sig2*gamma*U(1,1);

    for jj = 2:n
      beta = alpha;
      alpha = alpha + b(jj)*a(jj);
      lambda = -a(jj)*gamma;
      gamma = 1/alpha;
      U(jj,jj) = beta*gamma*U(jj,jj);
      for ii = 1:(jj-1)
        beta = U(ii,jj);
        U(ii,jj) = beta + b(ii)*lambda;
        b(ii) = b(ii) + b(jj)*beta;
      end
    end
    phat = phat0 + b'*r*gamma;
    SigmaP = mth_ud2p(U);
  else
    fprintf('\nInvalid update algorithm\n');
    phat = phat0;
    SigmaP = P0;
  end
