classdef GaussHermiteCubature < handle
    % Gauss-Hermite cubature sigma-points
    % 
    % DESCRIPTION
    %   Implements the multi-dimensional Gauss-Hermite cubature rules
    %   according to [1].
    %
    % PROPERTIES
    %   p (r/w, default: 3)
    %       Order of the polynomial to use.
    %
    % METHODS
    %   GaussHermiteCubature(p)
    %       Initializes the cubature. 'p' is the order of the polynomial to
    %       be used as described above.
    %
    %   [X, wm, wc] = calculateSigmaPoints(m, P)
    %       Calculates the sigma points using the mean m and covaraince P
    %       and returns the NxP^N matrix of sigma points (each column is a
    %       sigma point) together with the weights 'wm' for the mean and
    %       'wc' for the covariance. Note that wm and wc are the same for
    %       Gauss-Hermite cubatures.
    %
    % REFERENCES
    %   [1] S. S?rkk?, "Bayesian Filtering and Smoothing", Cambridge
    %       University Press, 2013
    %
    % SEE ALSO
    %   UnscentedTransform
    %
    % VERSION
    %   2016-12-19
    % 
    % AUTHORS
    %   Roland Hostettler <roland.hostettler@aalto.fi>   
    
    %#ok<*PROPLC>
    %#ok<*PROP>
        
    %% Properties
    properties (Access = public)
        % Order of the approximation
        p = 3;
    end
    
    %% Private Properties
    properties (Access = private)
        % Stores the roots and their 1D-weights
        xi;
        w;
        
        % Stores the unit sigma-points and their weights
        Xi;
        W;
    end
    
    %% Property Access Methods
    methods
        %% Setter for p
        % Also recalculates the roots of the polynomial whenever p is
        % altered
        function set.p(self, p)
            if p < 1
                error('Too low polynomial order.');
            end
            self.p = p;
            self.calculateRoots();
        end
    end

    %% Methods
    methods (Access = public)
        %% Constructor
        function self = GaussHermiteCubature(p)
            if nargin == 1
                self.p = p;
            end
            self.calculateRoots();
        end
        
        %% Calculate the Sigma Points
        function [X, wm, wc] = calculateSigmaPoints(self, m, P)
            Nx = size(m, 1);
            if Nx ~= size(self.Xi, 1)
                self.calculateUnitSigmaPoints(Nx);
            end
            wm = self.W;
            wc = self.W;
            X = m*ones(1, size(self.Xi, 2)) + chol(P).'*self.Xi;
        end
    end
    
    %% Private Methods
    methods (Access = private)
        %% Calculates the Unit Sigma-Points
        function calculateUnitSigmaPoints(self, N)
            xi = self.xi;
            w = self.w;
            p = self.p;
            Xi = zeros(N, p^N);
            W = zeros(1, p^N);

            % The goal is to create the Nxp^N matrix of sigma points Xi and
            % the 1xp^N weight vector. The matrix has the following form:
            %        _                                                                 _
            %       | xi(1) xi(1) ... xi(1) xi(2) xi(2) ... xi(2) xi(p) xi(p) ... xi(p) |
            %       | xi(1) xi(1)     xi(p) xi(1) xi(1) ... xi(p) xi(1) xi(1) ... xi(p) |
            %  Xi = | ...   ...       ...   ...   ...       ...   ...   ...      ...    |
            %       |_xi(1) xi(2) ... xi(p) xi(1) xi(2) ... xi(p) xi(1) xi(2) ... xi(p)_|
            %
            % The weights are constructed in the same way and in the end
            % the product of all the weights is taken.
            
            % Some debugging code below.
            % xi = [1, 2].';
            % w = [0.4, 0.6].';
            % p = length(xi);
            % N = 3;
            % Xi = zeros(N, p^N);
            
            for n = 1:N
                % First we create the basic repetitive pattern, i.e. how
                % many times x(1) should be repeated in this row.
                % Then, we repeat that pattern p^(n-1) times to build a
                % full row and put everything together.
                rxi = ones(p^(N-n), 1)*xi(:).';                
                rxi = rxi(:)*ones(1, p^(n-1));
                Xi(n, :) = rxi(:).';
                
                % Calculate the weights
                rw = ones(p^(N-n), 1)*w(:).';
                rw = rw(:)*ones(1, p^(n-1));
                W(n, :) = rw(:).';
            end
            W = prod(W, 1);
            
            self.Xi = Xi;
            self.W = W;
        end
        
        %% Gauss-Hermite Polynomial
        function y = calculateHermitePolynomial(self, p, x)
            switch p
                case 0
                    y = 1;
                case 1
                    y = x;
                case 2
                    y = x.^2 - 1;
                case 3
                    y = x.^3 - 3*x;
                case 4
                    y = x.^4 - 6*x.^2 + 3;
                otherwise
                    y = ( ...
                        x.*self.calculateHermitePolynomial(p-1, x) ...
                        - (p-1)*self.calculateHermitePolynomial(p-2, x) ...
                    );
            end
        end
        
        %% Calculates the Roots of the Gauss-Hermite Polynomial of Order p
        function calculateRoots(self)
            p = self.p;
            switch p
                case 0
                    xi = [];
                case 1
                    xi = 0;
                case 2
                    xi = [-1, 1].';
                case 3
                    xi = [-sqrt(3), 0, sqrt(3)].';
                case 4
                    xi = [
                        -sqrt((6+sqrt(24))/2);
                        -sqrt((6-sqrt(24))/2);
                         sqrt((6-sqrt(24))/2);
                         sqrt((6+sqrt(24))/2);
                    ];
                otherwise
                    error('Not implemented yet.');
            end
            self.w = factorial(p)./(p^2*self.calculateHermitePolynomial(p-1, xi).^2);
            self.xi = xi;
        end
    end
end
