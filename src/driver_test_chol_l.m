
clear;

P = [ 3 2 1 ; 2 4 3 ; 1 3 5 ];

S = mth_sqrtm(P);
[~, R] = mth_qr(S');

fprintf('\nP - SST:  %1.3e', max(max(abs(P - S*S'))));
fprintf('\nP - RTR:  %1.3e', max(max(abs(P - R'*R))));

c = .25;
x = [.1 .3 .2]';
P2 = P + c*x*x';

L = R';
L2 = mth_chol_upd_l(L, c, x);
fprintf('\nUpdate');
fprintf('\nP - LLT:  %1.3e', max(max(abs(P2 - L2*L2'))));

P2 = P - c*x*x';
L2 = mth_chol_upd_l(L, -c, x);
fprintf('\nDowndate');
fprintf('\nP - LLT:  %1.3e', max(max(abs(P2 - L2*L2'))));

fprintf('\n');

