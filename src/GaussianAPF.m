classdef GaussianAPF < handle
    % Gaussian approximation auxiliary particle filter
    % 
    % DESCRIPTION
    %   
    %
    % PROPERTIES
    % 
    %
    % METHODS
    %
    %
    % SEE ALSO
    %
    %
    % VERSION
    %   2017-01-24
    % 
    % AUTHORS
    %   Roland Hostettler <roland.hostettler@aalto.fi>   
    
    % TODO
    %   * Implement more generic moment matching
    %   * Implement sanity checks
    %   * Implement handling of control input
    %   * Add complete documentation
    %   * Move to yamfast
    %   * (Maybe) Move calculateMoments out of class
    
    % Disable some annoying warnings
    %#ok<*PROPLC>
    
    %% Public Properties
    properties (Access = public)
        % No. of particles, default to 100
        M = 100;
        
        % Particles & weights
        x;
        w;
        
        % Linearization approach for moment matching
        linearization = 'ekf';
        
        % No. of posterior linearization iterations
        J = 0;
        
        % Estimated state
        xhat;
        
        % A Wiener State-Space Model
        model;
    end
    
    %% Private Properties
    properties (Access = private)
        % Sigma-point rule
        rule = UnscentedTransform();
    end

    %% Methods
    methods (Access = public)
        %% Initialize the Filter
        function self = GaussianAPF(model, M, type)
            % Store the model
            self.M = M;
            self.model = model;
            
            % Initialize
            Nx = size(model.m0, 1);
            self.x = model.m0*ones(1, M) ...
                + chol(model.P0, 'lower')*randn(Nx, M);
            self.w = 1/M*ones(1, M);
            
            % Linearization method
            if nargin == 3
                if ischar(type) && strcmp(type, 'ekf')
                    self.linearization = 'ekf';
                elseif ~isempty(type)
                    self.linearization = 'sigma-point';
                    self.rule = type;
                else
                    error('Something is wrong with the chosen linearization.');
                end
            end
        end
        
        %% Filter Iteration
        function xhat = update(self, y, t, u)
            %% Preparations & quick access
            x = self.x;
            w = self.w;
            [Nx, M] = size(x);
            Ny = size(y, 1);

            % Preallocate
            xp = zeros(Nx, M);
            Pp = zeros(Nx, Nx, M);
            xn = zeros(Nx, M);
            yp = zeros(Ny, M);
            B = zeros(Nx, Ny, M);
            S = zeros(Ny, Ny, M);
            mu_x = zeros(Nx, M);
            Sigma_x = zeros(Nx, Nx, M);
    
            %% Calculate the Importance Distribution's moments
            for m = 1:M
                [xp(:, m), Pp(:, :, m), yp(:, m), S(:, :, m), B(:, :, m)] = self.calculateProposalMoments(y, x(:, m), t, u);
                mu_x(:, m) = xp(:, m) + B(:, :, m)/S(:, :, m)*(y - yp(:, m));
                Sigma_x(:, :, m) = Pp(:, :, m) - B(:, :, m)/S(:, :, m)*B(:, :, m)';
            end
            
            %% Resample
            % Auxiliary variable probabilities
            % TODO: Implement this in log-space
            v = w.*mvnpdf(y.', yp.', S).';
            v = v/sum(v);
            alpha = resample(v);
            
            %% Draw New Particles
            for m = 1:M
                xn(:, m) = mu_x(:, alpha(m)) ...
                    + chol(Sigma_x(:, :, alpha(m)), 'lower')*randn(Nx, 1);
                
                % TODO: Implement in log-space
                mu = yp(:, alpha(m)) ...
                    + B(:, :, alpha(m))'/Pp(:, :, alpha(m))*(xn(:, m) - xp(:, alpha(m)));
                Sigma = S(:, :, alpha(m)) ...
                    - B(:, :, alpha(m))'/Pp(:, :, alpha(m))*B(:, :, alpha(m));
                w(:, m) = self.model.py_eval(y, xn(:, m), t, u)/mvnpdf(y, mu, Sigma).';
            end
            w = w/sum(w);
        
            %% Estimate & store results
            % TODO: store covariance, map
            xhat = xn*w';
            self.xhat = xhat;
            self.x = xn;
            self.w = w;
        end
        
        %% Filter a Batch of Data
        function xhat = filter(self, y, t, u)
            N = size(y, 2);
            if nargin <= 2 || isempty(t)
                t = 1:N;
            end
            if nargin <= 3 || isempty(u)
                u = zeros(1, N);
            end
            xhat = zeros(size(self.x, 1), N);
            
            for n = 1:N
                xhat(:, n) = self.update(y(:, n), t(n), u(n));
                self.model.t = t(n);
            end
        end
    end
    
    %% Helper Methods
    methods (Access = protected)
        %% Moment Matching
        % Calculates the moments of the Gaussian approximation
        %                                   _      _    _    _   _     _
        %       p(x[n], y[n] | x[n-1]) = N(|  x[n]  |  |  xp  | | Q  D  |)
        %                                  |_ y[n] _|, |_ yp _| | D' S _|
        function [xp, Q, yp, S, D] = calculateProposalMoments(self, y, x, t, u)
            model = self.model;
            
            % Prior Linearization
            Q = model.Q(x, t, u);
            xp = model.f(x, zeros(size(Q, 1), 1), t, u);
            [yp, S, D] = self.calculateSLRMoments(xp, Q, t, u);
            
            % Posterior linearization; set J to 0 for none (default)
            for j = 1:self.J
                % Calculate posterior
                K = D/S;
                mu_pi = xp + K*(y - yp);
                Sigma_pi = Q - K*S*K';

                % Calculate SLR moments
                [ybar, Sbar, Dbar] = self.calculateSLRMoments(mu_pi, Sigma_pi, t, u);

                % Calculate linear approximation
                An = Dbar/Sbar;
                bn = ybar - An*mu_pi;
                Sigma_v = Sbar - An*Sigma_pi*An';

                % Calculate updated moments
                yp = An*xp + bn;
                S = An*Q*An' + Sigma_v;
                D = Q*An';
            end
        end
        
        %% Calculates the oments used for SLR 
        function [yp, S, D] =  calculateSLRMoments(self, mu, Sigma, t, u)
            model = self.model;
            
            switch self.linearization
                case 'ekf'
                    R = model.R(mu, t, u);
                    [yp, Gx, Gr] = model.g(mu, zeros(size(R, 1), 1), t, u);
                    S = Gx*Sigma*Gx' + Gr*R*Gr';
                    D = Sigma*Gx';
                    
                case 'sigma-point'
                    % TODO: Assumes additive measurement noise for now
                    % Determine the state, noise, and measurement
                    % dimensions
                    Nx = size(mu, 1);
                    R = model.R(mu, t, u);
                    Nr = size(R, 1);
                    yp = model.g(mu, zeros(Nr, 1), t, u);
                    Ny = size(yp, 1);

                    [X, wm, wc] = self.rule.calculateSigmaPoints(mu, Sigma);
                    L = size(X, 2);
                    Y = zeros(Ny, L);
                    for j = 1:L
                        Y(:, j) = model.g(X(:, j), zeros(Nr, 1), t, u);
                    end

                    % Calculate the measurement and covariance
                    yp = Y*wm';
                    S = zeros(Ny);
                    D = zeros(Nx, Ny);
                    for j = 1:L
                        S = S + wc(j)*(Y(:, j) - yp)*(Y(:, j) - yp)';
                        D = D + wc(j)*(X(1:Nx, j) - mu)*(Y(:, j) - yp)';
                    end
                    S = S + R;
                                        
                otherwise
                    error('Unknown quadrature method: %s.', self.linearization);
            end
        end
    end
end
