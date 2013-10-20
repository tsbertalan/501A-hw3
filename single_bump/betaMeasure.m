function beta = betaMeasure(u, S)
    N = numel(u);
    indices = 1:N;
    o = ones(N, 1);
    P = permMat(S);
    v = P * reshape(u, 1, N)';
    v = v - min(v);
    s = N/10;
    
    top = max(v);
    bot = min(v);
    range = top - bot;
    
    center = findCenter(v);
    [closest, centerInd] = min(abs(v - center));
    inlocs = (centerInd-s < indices) & (indices < centerInd+s);

    totMass = sum(v);
    inMass = sum(v(inlocs));
    
    beta = 1 - (inMass / totMass);
    
    % % Demo
%     plot([centerInd, centerInd], [-1, 1], 'k')
%     plot([centerInd-s, centerInd-s], [-1, 1], 'r')
%     plot([centerInd+s, centerInd+s], [-1, 1], 'r')
end