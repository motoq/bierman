clear;

P = [ 3 2 1 ; 2 4 3 ; 1 3 5 ];

S = mth_sqrtm(P);
Sc = mth_chol(P);
[~, R] = mth_qr(S');                             % S = mth_qr(S)

fprintf('\nP - SST:  %1.3e', max(max(abs(P - S*S'))));
fprintf('\nP - ScTSc:  %1.3e', max(max(abs(P - Sc'*Sc))));

c = .25;
x = [.1 .3 .2]';
P2 = P + c*x*x';
R2 = cholupdate(R, sqrt(c)*x, '+');

fprintf('\nP2 - R2TR2:  %1.3e', max(max(abs(P2 - R2'*R2))));

S2 = mth_chol_upd(S, c, x);
fprintf('\nP2 - S2S2T:  %1.3e', max(max(abs(P2 - S2*S2'))));


P3 = P + P2;
S = mth_sqrtm(P);
S2 = mth_sqrtm(P2);
A = [S S2];
[~, R] = mth_qr(A');
fprintf('\nP3 - RqrTRqr:  %1.3e', max(max(abs(P3 - R'*R))));

R = mth_householder_tri(A');
fprintf('\nP3 - RhTRh:  %1.3e', max(max(abs(P3 - R'*R))));


A = [S ; S2];

[~, R] = mth_qr(A);
R = R(1:3,1:3);
fprintf('\nP3 - Rqr2Rqr2T:  %1.3e', max(max(abs(P3 - R*R'))));

R = mth_householder_tri(A);
R = R(1:3,1:3);
fprintf('\nP3 - Rh2Rh2T:  %1.3e', max(max(abs(P3 - R*R'))));

fprintf('\n');