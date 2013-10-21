clear all;
close all;
figure();
hold all;
styles = {'-k', '--k'};

Ncells = 250;
indices = 1:1:Ncells;
sigma_noise = 0.4;

dW = .1;
wDlist = -.5:dW:.5;
rMaxList = 1:1:numel(wDlist);

for s=[0, 1]
    initialU = s * 4 * exp(-(indices - Ncells / 2).^2 / ...
        (2 * sigma_noise.^2));

    for i=1:1:numel(wDlist) 
        disp(wDlist(i))
        [r, u] = single_bump('T', 20, 'wDA', wDlist(i), 'do_plot', 1, ...
                             'Ncells', Ncells, 'initialU', initialU, ...
                             'sigma_noise', sigma_noise, 'use_field_B', 0);
        rMaxList(i) = max(r);
    end
    plot(wDlist, rMaxList, char(styles(s+1)));
    xlabel('w_D');
    ylabel('max(r)');
end
legend({'without initial bump', 'with initial bump'});
saveas(gcf(), sprintf('initialbump-dW_%.2f.eps', dW));