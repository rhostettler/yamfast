function px = logmvnpdf(x, m, C)
% Logarithm of Multivariate Normal PDF
%
% SYNOPSIS
%   px = logmvnpdf(x)
%   px = logmvnpdf(x, m, C)
%
% DESCRIPTION
%   Returns the logarithm of N(x; m, C) or N(x; 0, I) if m and C are
%   omitted, i.e. the log-likelihood. Everything is calculated in
%   log-domain such that numerical precision is retained.
%
%   The arguments x and m are automatically expanded to match each other.
%
% PARAMETERS
%   x   MxN vector of values to evaluate.
%   m   MxN vector of means (optional, default: 0).
%   C   Covariance matrix (optional, default: I).
%
% VERSION
%   2016-12-02
%
% AUTHORS
%   Roland Hostettler <roland.hostettler@aalto.fi>

% TODO:
%   * Include sanity checks for C and the M-sizes of the vectors.
%   * Consider implementation using Cholesky-factor instead.

    %% Autocomplete
    switch nargin
        case 1
            m = zeros(size(x));
            C = eye(size(x, 2));
        case 2
            C = eye(size(x, 2));
    end

    % Some sanity checks
    [Nx, Mx] = size(x);
    [Nm, Mm] = size(m);
    Nv = size(C, 1);
    
    % Automagically expand the arguments
    if Nx ~= Nm
        if Nx == 1 && Nm > 1
            x = ones(Nm, 1)*x;
        elseif Nx > 1 && Nm == 1
            m = ones(Nx, 1)*m;
        else
            error('Input argument size mismatch');
        end
    end
    
    %% Calculation
%     Cinv = C\eye(Nv);
%     px = -Nv/2*log(2*pi) - 1/2*log(det(C)) - 1/2*diag((x-m)*Cinv*(x-m)');
    L = chol(C);
    epsilon = (x-m)/L;
    px = -Nv/2*log(2*pi) - 1/2*log(det(C)) - 1/2*diag(epsilon*epsilon');
end
