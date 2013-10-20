function P = permMat(B_order)
    P = eye(numel(B_order));
    P = P(:,B_order);  % permutation matrix for the A -> B map
end