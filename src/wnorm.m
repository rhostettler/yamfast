function n = wnorm(x, W)
% Calculates the weighted norm x'*W*x
%
% SYNOPSIS
%   n = wnorm(x, P)
%
% DESCRIPTION
%   Calculates the weighted norm x^T W x in a numerically stable way.
%   Additionally, Nx x M input arrays for x are accepted, calculating the
%   norm of multiple vectors in one go.
%
% PARAMETERS
%   x   NxM matrix where each column corresponds a vector.
%   W   NxN weight matrix.
%
% RETURNS
%   n   1xM vector of weighted norms.
%
% VERSION
%   2017-01-17
%
% AUTHORS
%   Roland Hostettler <roland.hostettler@aalto.fi>

    %% Sanity Checks
    switch nargin
        case 1
            W = eye(size(x, 1));
        case 2
            % nop
        otherwise
            error('Error using wnorm, see help wnorm.');
    end

    %% Calculations
    L = chol(W);        % W = L'*L
    f = L*x;
    n = diag(f'*f).';   % (L*x)'*(L*x) = x'*L'*L*x = x*W*x
end
