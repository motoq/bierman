function [x, y, z] = mth_helix_parametric(t, r, c)
% MTH_HELIX_PARAMETRIC generates a helix as a function of an independent
% variable.  The helix is centered about the z-axis.
%
%-----------------------------------------------------------------------
% Copyright 2018 Kurt Motekew
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%-----------------------------------------------------------------------
%
% Inputs:
%   t   Independent parameter, range of angles over which the curve is to
%       be generated, radians
%   r   Radius of helix
%   c   Rate at which to extend the helix along the z-axis   
%
% Return:
%   x  Cartesian x coordinates, in units of r
%   y  Cartesian y coordinates, in units of r
%   z  Cartesian z coordinates, in units of r
%
% Kurt Motekew   2018/07/19
%

  x = r*cos(t);
  y = r*sin(t);
  z = c*t;

