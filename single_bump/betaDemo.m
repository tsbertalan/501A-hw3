% betaMeasure demo
close all; 
N = 50;

noises = [0, .5]
for i=1:2
    figure()
hold all;

    noise = noises(i)
    loc = locs(i)
    
u = makeBump(N, loc, 1, 2) + rand(1,N)*noise;
s = N/10;
plot(u);
5
center = findCenter(u);
[closest, centerInd] = min(abs(u - center));

verti(centerInd, 'k');
verti(centerInd-s, 'r');
verti(centerInd+s, 'r');

xlabel('S')
ylabel('u_S')
title(sprintf('\\beta_S=%.2f', betaMeasure(u, 1:N)))
saveas(gcf(), sprintf('betaMeasure-noise_%.2f.eps', noise), 'epsc');
end

% legend({...
%     sprintf('\\mu_0=%d', 8), ...
%     sprintf('\\mu_0=%d', N-8), ...
%     })

