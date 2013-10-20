function p2q7()
%P2Q6 Mixed stimulus

clear all; close all;

wDA = 0.3;
wDB = 0.3;

Ncells = 100;
dVA = -1;
dVB = 1;
gI = 0.3;
wE = 0.8;
[r, u] = single_bump('wDA', wDA, 'wDB', wDB, 'Ncells', Ncells,...
                     'dVA', dVA, 'dVB', dVB, 'gI', gI, 'wE', wE)


end
