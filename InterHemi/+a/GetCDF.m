function [cdf, bin] = GetCDF(vec, binsize)
%
% creates a cdf by cumsuming a hist of vec with binsize
%
% arb 040215
if isnan(vec)
    cdf = 0;
    bin = 0;
    return;
end
[count, bin] = hist(vec, range(vec)/binsize);

cdf = cumsum(count)/sum(count);