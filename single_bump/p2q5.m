function [ output_args ] = p2q5( input_args )
%P2Q5 Mixed initial condidions.
%   A mix between belief in space A and space B.
clear all;
close all;
Ncells = 50;
sigmaE = 4;
B_order = randperm(Ncells);
figure()
hold all;

initialA = makeBump(Ncells, 5, 1, sigmaE) * .8;
plot(initialA)
initialB = makeBump(Ncells, 45, 1, sigmaE) * .1;
plot(initialB)
P = permMat(B_order);
initialB = (P * initialB')';
initialU = initialA + initialB;
plot(initialU)


single_bump('initialU', initialU, 'B_order', B_order, 'do_plot', 1, ...
            'Ncells', Ncells)


end

