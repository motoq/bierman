clear;

P = [ 3 2 1 ; 2 4 3 ; 1 3 5 ];

S = mth_sqrtm(P);
[~, R] = mth_qr(S');                             % S = mth_qr(S)

c = .25;
x = [.1 .3 .2]';
P2 = P + .25*x*x';
R2 = cholupdate(R, sqrt(c)*x, '+');

fprintf('\nP2 - R2TR2:  %1.3e\n', max(max(abs(P2 - R2'*R2))));

S2 = mth_chol_upd(S, c, x);
fprintf('\nP2 - W2S2T:  %1.3e\n', max(max(abs(P2 - S2*S2'))));