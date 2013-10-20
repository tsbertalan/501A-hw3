function u = makeBump( Ncells, u0, scale, sigmaE)
%MAKEBUMP make an activity bump
indices = 1:Ncells;
u = scale * exp(-(indices - u0).^2 / ...
            (2 * sigmaE.^2));
end

