function [ center ] = findCenter( rho_ )
    %FINDCENTER locate the center-of-mass of a sequence of values
    %   A weighted mean of indices. Provided as a float; round to get
    %   a corresponding index into the original data.
    lowest = min(rho_);
    N = numel(rho_);
    rho = reshape(rho_ - lowest, N, 1);
    M = sum(rho);
    r = 1:1:numel(rho);
    center = sum(r * rho) / M;
%     center = round(center);
end