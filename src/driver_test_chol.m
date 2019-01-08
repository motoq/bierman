clear;

P = [ 3 2 1 ; 2 4 3 ; 1 3 5 ];

S = mth_sqrtm(P);
[~, R] = mth_qr(S');                             % S = mth_qr(S)

c = .25;
x = [.1 .3 .2]';
P2 = P + c*x*x';
R2 = cholupdate(R, sqrt(c)*x, '+');

fprintf('\nP2 - R2TR2:  %1.3e\n', max(max(abs(P2 - R2'*R2))));

S2 = mth_chol_upd(S, c, x);
fprintf('\nP2 - S2S2T:  %1.3e\n', max(max(abs(P2 - S2*S2'))));


P3 = P + P2;
A = [S S2];
[~, R3] = mth_qr(A');
fprintf('\nP3 - R3TR3:  %1.3e\n', max(max(abs(P3 - R3'*R3))));

Rh = mth_householder_tri(A');
fprintf('\nP3 - RhTRh:  %1.3e\n', max(max(abs(P3 - Rh'*Rh))));

