function [vel, SigmaV] = traj_pos2vel(dt, pos1, SigP1, pos2, SigP2,...
                                                        pos3, SigP3)
% TRAJ_POS2VEL Estimates velocity given three position estimates, associated
% uncertainties, and the time interval between each position.
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
%   dt      Time between position values
%   pos1    Initial position, [3x1]
%   SigP1   Initial position covariance, [3x3]
%   pos1    Middle position, [3x1]
%   SigP1   Middle position covariance, [3x3]
%   pos3    Final position, [3x1]
%   SigP3   Final position covariance, [3x3]
%
% Return:
%   vel      Estimated velocity at current time, [3x1]
%   SigmaV   Estimated velocity covariance, [3x3]
%
% Kurt Motekew   2014/11/19
%

W1 = (SigP1 + SigP2)^-1;
W2 = (SigP2 + SigP3)^-1;
Ap = dt*eye(3);

SigmaV = (Ap'*W1*Ap + Ap'*W2*Ap)^-1;
vel = SigmaV*(Ap'*W1*(pos2-pos1) + Ap'*W2*(pos3-pos2));
