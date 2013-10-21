clear all;
close all;

Ncells = 50;
indices = 1:1:Ncells;

du0 = .1;
width = 3;
T = 50;
u0list = Ncells/2-width:.1:Ncells/2+width;

wD = 0;
sigmaE = 4;

styles = {'-k', '-k'};
noises = [0   , .4];  % Dictionaries would be nice about now.
f = figure('paperposition', [0, 0, 4.0, 3.0]);
hold all;

% % for animation
% g = figure()

plots = [];
legends = {};

for i=1:1:numel(noises)
    sigma_noise = noises(i);
    style = styles{i};
    movements = 1:1:numel(u0list);
    for j=1:1:numel(u0list)
        u0 = u0list(j);
        
        initialU = makeBump(Ncells, u0, 4, sigmaE);

        [r, u] = single_bump('T', T, 'wDA', wD, 'do_plot', 0, 'use_field_B', 0,...
                             'Ncells', Ncells, 'initialU', initialU, ...
                             'sigma_noise', sigma_noise, 'sigmaE', sigmaE);
        uf = findCenter(u);
     
        % % pseudo-animation!        
%         figure(g)
%         clf
%         hold all;
%         plot(initialU, '-k')
%         plot(u, '--k')
%         plot([u0, u0], [2, 4], '-k');
%         plot([uf, uf], [-2, 2], '--k');
%         legend({'initial', 'final'})

        movements(j) = uf - u0;

    end
    figure()
    hold all;
    handles(i) = plot(u0list, movements, style);
    title(sprintf('\\sigma_\\eta=%.1f', sigma_noise))
    xlabel('u(0)')
    ylabel(sprintf('u(%d) - u(0)', T))
    horiz(0, 'k')
    
    saveas(gcf(), sprintf('movement-sn_%.2f.eps', sigma_noise), 'epsc');
end



