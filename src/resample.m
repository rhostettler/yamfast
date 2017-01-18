function i = resample(w)
% Systematic resampling
%
% SYNOPSIS
%   i = resample(w)
%
% DESCRIPTION
%   Systematic resampling, returns randomized indices i such that
%   Pr(i) = w(i).
%
% PARAMETERS
%   w   Probabilities.
%
% RETURNS
%   i   The indices.
%
% VERSION
%   2017-01-17
%
% AUTHORS
%   Roland Hostettler <roland.hostettler@aalto.fi>

    w = w/sum(w);
    M = length(w);
    i = zeros(1, M);
    k = 0;
    u = 1/M*rand();
    for j = 1:M
        N = floor(M*(w(j)-u)) + 1;
        if N > 0
            i(k+1:k+N) = j;
            k = k + N;
        end
        u = u + N/M - w(j);
    end
    i = i(randperm(M));
end
