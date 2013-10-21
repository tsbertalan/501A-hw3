function [u, r] = p2q5(  )
%P2Q5 Mixed initial condidions.
%   A mix between belief in space A and space B.
clear all;
close all;
Ncells = 50;
sigmaE = 4;
T = 30;
trials = 1024;

Afrac = 1.0;
Bfrac = 1.0;

uA = 32;
uB = 1;

betaPoints = zeros(trials, 2);

initialA = makeBump(Ncells, uA, 1, sigmaE) * Afrac;
initialB = makeBump(Ncells, uB, 1, sigmaE) * Bfrac;


for trial=1:trials
    disp(trial);
    B_order = randperm(Ncells);
    P = permMat(B_order);
    initialB = (P * initialB')';
    initialU = initialA + initialB;
    
    [r, u] = single_bump('initialU', initialU, 'B_order', B_order, ...
                         'do_plot', 0, 'Ncells', Ncells, 'T', T);
    
    bA = betaMeasure(r, 1:Ncells);
    bB = betaMeasure(r, B_order);
    
    % find cases where either peak is strong
    thresh = .12;
    if bA < thresh & bB < thresh
        figure(19);
        clf;
        hold all;
        plot(r, 'r');
        plot(P * r, 'b');
        xlabel('indexing in A or B space');
        ylabel('rate in A or B space');
        legend({'u_A', 'u_B'});
        title(sprintf('\\beta_A=%.1f, \\beta_B=%.1f',...
                      [betaMeasure(r, 1:Ncells), betaMeasure(r, B_order)]))
        verti(findCenter(r), '--r');
        verti(findCenter(P * r), '--b');
        disp(sprintf('bA=%.2f ; bB=%.2f', [bA, bB]));
        
        if bA < thresh & bB < thresh
            saveas(gcf(), sprintf('p2q5-strongPeaks-bA_%.2f-bB_%.2f.eps', [bA, bB]), 'epsc');
        end

    end
    betaPoints(trial, :) = [bA, bB];
end

figure();
hold all;
scatter(betaPoints(:,1), betaPoints(:,2));
xlim([0, 1]);
ylim([0, 1]);
xlabel('\beta_A');
ylabel('\beta_B');
title(sprintf('\\mu_A=%.1f, \\mu_B=%.1f', [uA, uB]));

m = mean(betaPoints, 1);
bA = m(1);
bB = m(2);
verti(bA, 'k');
horiz(bB, 'k');

saveas(gcf(), sprintf('p2q5-betaScatter-muA_%.1f-muB_%.1f-%d_trials', [uA, uB, trials]), 'epsc')
end