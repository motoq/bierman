function matrix3X3_plot(A, N, makesurface)
%
% Plots a 3x3 covariance matrix.
%
% Inputs:
%   A            Covariace matrix to plot - 3x3
%   N            Number of something - see 'help ellipsoid' for its
%                useless documentation
%   makesurface  If true, create a surface.  Otherwise, a mesh
%
% Kurt Motekew   2014/10/19
%
 
[XX, YY, ZZ] = matrix3X3_points(A,N);

if (makesurface)
  surf(XX,YY,ZZ);
else
  mesh(XX,YY,ZZ);
  colormap([.1 .2 .3]);
  hidden('off');
end
xlabel('x');
ylabel('y');
zlabel('z');
axis equal;
 
end
