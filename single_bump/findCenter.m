function [ center ] = findCenter( rho_ )
    %FINDCENTER locate the center-of-mass of a sequence of values
    %   A weighted mean of indices. Provided as a float; round to get
    %   a corresponding index into the original data.
    lowest = min(rho_);
    rho = rho_ - lowest;
    M = sum(rho);
    r = 1:1:numel(rho);
    center = sum(r * rho) / M;
%     center = round(center);
end