function p2q6()
%P2Q6 noiseless double-enviroment network.

% close all;
% clear all;

Ncells = 50;
T = 20;
sigmaE = 4;
sigma_noise = 0;

u0s = 1:Ncells;

movements = [];
B_order = randperm(Ncells);

for i=1:Ncells
    u0 = u0s(i)
    initialU = makeBump(Ncells, u0, 1, sigmaE);
    [r, u] = single_bump('initialU', initialU, 'T', T, 'Ncells', Ncells,...
                         'sigmaE', sigmaE, 'sigma_noise', sigma_noise,...
                         'do_plot', 0, 'B_order', B_order);
     
     center = findCenter(r);
     movements(i) = u0 - center;
end

figure()
hold all;
scatter(u0s, movements);
xlabel('u_0');
ylabel('u_f - u_0');
horiz(0, 'k');
saveas(gcf(), 'p2q6.eps', 'epsc');

end

