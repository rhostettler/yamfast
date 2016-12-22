classdef UnscentedTransform < handle
    % Unscented transform (UT) (and cubature) sigma-points
    % 
    % DESCRIPTION
    %   Implements the unscented transform (UT) [1] as it is used in, for
    %   example, the unscented Kalman filter or approximations of the
    %   optimal proposal in particle filters.
    %
    %   Note that by setting alpha = 1, beta = 0, kappa = 0, the
    %   cubature integration rule according to [2] is obtained.
    %
    % PROPERTIES
    %   alpha (r/w, default: 1)
    %       'alpha'-parameter of the UT.
    %
    %   beta (r/w, default: 0)
    %       'beta'-parameter of the UT.
    %
    %   kappa (r/w, default: 1)
    %       'kappa'-parameter of the UT.
    %
    % METHODS
    %   UnscentedTransform(alpha, beta, kappa)
    %       Initializes the UT. The parameters alpha, beta, and kappa are
    %       the user-specified parameters as described above. These are
    %       optional and if omitted, the default values above are used.
    %
    %   [X, wm, wc] = calculateSigmaPoints(mu, Sigma)
    %       Calculates the sigma points based on the mean 'mu' and
    %       covariance matrix 'Sigma' and returns the sigma-points in 'X'
    %       where each column corresponds to one sigma-point and with
    %       weights 'wm' (mean) and 'wc' (covariance).
    %
    % REFERENCES
    %   [1] Wan, E. A. and Van der Merwe, R. "The unscented Kalman filter", 
    %       in Kalman Filtering and Neural Networks, Wiley, 2001
    %
    %   [2] Arasaratnam, I. and Haykin, S. "Cubature Kalman filters", IEEE 
    %       Transactions on Automatic Control, 54(6), 1254?1269, 2009
    %
    % SEE ALSO
    %   GaussHermiteCubature
    %
    % VERSION
    %   2016-12-18
    % 
    % AUTHORS
    %   Roland Hostettler <roland.hostettler@aalto.fi>   
    
    %% Properties
    properties (Access = public)
        % Spread
        alpha = 1;
        
        % Secondary spread
        kappa = 1
        
        % Prior knowledge
        beta = 0
    end
    
    %% Methods
    methods (Access = public)
        %% Constructor
        function self = UnscentedTransform(alpha, beta, kappa)
            if nargin >= 1 && ~isempty(alpha)
                self.alpha = alpha;
            end
            if nargin >= 2 && ~isempty(beta)
                self.beta = beta;
            end
            if nargin == 3 && ~isempty(kappa)
                self.kappa = kappa;
            end
        end
        
        %% Calculates the Sigma-Points
        function [X, wm, wc] = calculateSigmaPoints(self, mu, Sigma)
            M = size(mu, 1);
            lambda = self.alpha^2*(M + self.kappa) - M;
            wm = 1/(2*(M+lambda))*ones(1, 2*M+1);
            wm(1) = lambda/(M+lambda);
            wc = wm;
            wc(1) = lambda/(M+lambda) + (1-self.alpha^2+self.beta);
            dx = sqrt(M+lambda)*chol(Sigma, 'lower');
            X = [mu, mu*ones(1, M) + dx, mu*ones(1, M) - dx];
        end
    end
end
