function Phi = traj_strans(dt)
% TRAJ_STRANS Returns the state transition matrix for a simple flat earth
% ballistic trajectory where position and velocity form the state vector.
% This form omits any contribution due to gravitational acceleration.
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
%   dt   Change in time from 
%
% Return:
%   Phi   [3x3] State transition matrix, dx_new/dx_old, x = [pos ; vel]
%
% Kurt Motekew   2016/11/22
%

  Phi = eye(6);
  Phi(1,4) = dt;
  Phi(2,5) = dt;
  Phi(3,6) = dt;
