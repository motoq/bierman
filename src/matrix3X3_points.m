function [XX, YY, ZZ] = matrix3X3_points(A,N)
%
% Creates x, y, and z points suitable for surf and mesh plots given a
% 3x3 covariance matrix.  This is based on a program from Nima Moshtagh
% with some corrections.
%
% Inputs:
%   A            Covariace matrix to plot - 3x3
%   N            Number of something - see 'help ellipsoid' for its
%                useless documentation
% Return:
%   XX           X-axis points
%   YY           Y-axis points
%   ZZ           Z-axis points
%
%  Nima Moshtagh Original code, 2009, see below
%  Kurt Motekew  MODS
%    20140312  - Split from matrix_plot - only generates points - no plotting
%              - Correction such that axes are the square root of the
%                eigenvalues, not the eigenvalues.
%              - Changed doc a little and actually made it useful
%
%------------------------------------------------------------------------------
% Copyright (c) 2009, Nima Moshtagh
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
%
%    * Redistributions of source code must retain the above copyright 
%      notice, this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright 
%      notice, this list of conditions and the following disclaimer in 
%      the documentation and/or other materials provided with the distribution
%      
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
%AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
%IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
%ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
%LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
%CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
%SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
%INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
%CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
%ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
%POSSIBILITY OF SUCH DAMAGE.
%------------------------------------------------------------------------------
[~, D, V] = svd(A);
 
%----------------------------------
% generate the ellipsoid at (0,0,0)
%----------------------------------
a = sqrt(D(1,1));
b = sqrt(D(2,2));
c = sqrt(D(3,3));
 
[X,Y,Z] = ellipsoid(0,0,0,a,b,c,N);
    
%-------------------------
% rotate and the ellipsoid
%-------------------------
 
XX = zeros(N+1,N+1);
YY = zeros(N+1,N+1);
ZZ = zeros(N+1,N+1);
 
for k = 1:length(X),
  for j = 1:length(X),
    point = [X(k,j) Y(k,j) Z(k,j)]';
    P = V * point;
    XX(k,j) = P(1);
    YY(k,j) = P(2);
    ZZ(k,j) = P(3);
  end
end
