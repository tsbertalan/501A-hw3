function [u, r] = p2q5(  )
%P2Q5 Mixed initial condidions.
%   A mix between belief in space A and space B.
clear all;
close all;
Ncells = 50;
sigmaE = 4;
T = 30;
trials = 50;

Afrac = 1.0;
Bfrac = 1.0;

uA = 10;
uB = 40     ;

betaPoints = zeros(trials, 2);


figure(2);
hold all;

initialA = makeBump(Ncells, uA, 1, sigmaE) * Afrac;
plot(initialA)
initialB = makeBump(Ncells, uB, 1, sigmaE) * Bfrac;
plot(initialB)


for trial=1:trials
    disp(trial)
    B_order = randperm(Ncells);
    P = permMat(B_order);
    initialB = (P * initialB')';
    initialU = initialA + initialB;
    

%     figure(1)
%     hold all;
    [r, u] = single_bump('initialU', initialU, 'B_order', B_order, ...
                         'do_plot', 0, 'Ncells', Ncells, 'T', T);
%     plot(u, '-k');
    bA = betaMeasure(r, 1:Ncells);
    bB = betaMeasure(r, B_order);
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

figure();
hold all;
plot(r, 'r');
plot(P * r, 'b');
legend({'u_A', 'u_B'});
title(sprintf('\\beta_A=%.1f, \\beta_B=%.1f',...
              [betaMeasure(r, 1:Ncells), betaMeasure(r, B_order)]))
verti(findCenter(r), 'r');
verti(findCenter(P * r), 'b');

end