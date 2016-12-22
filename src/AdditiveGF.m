classdef AdditiveGF < handle
    % Additive Gaussian assumed density filter using sigma-points
    % 
    % DESCRIPTION
    %   Gaussian assumed density filter for non-linear state space models
    %   with additive process- and measurement noise of the form
    %
    %       x[n] = f(x[n-1], t[n], u[n]) + q[n]
    %       y[n] = g(x[n], t[n], u[n]) + r[n]
    %
    %   with x[0] ~ N(m[0], P[0]), q[n] ~ N(0, Q[n]), r[n] ~ N(0, R[n]).
    %
    %   The filter uses sigma-points to approximate the moment matching
    %   integrals. Typical filters of this class are:
    %
    %       * Unscented Kalman filter (UKF)
    %       * Cubature Kalman filter (CKF)
    %       * Gauss-Hermite Kalman filter (GHKF)
    %
    %   By default, the unscented transform is used, that is, the filter is
    %   an unscented Kalman filter.
    %
    % PROPERTIES
    %   TODO: Document these.
    %
    % METHODS
    %   TODO: Document these.
    %
    % REFERENCSE
    %   [1] S. Särkkä, "Bayesian Filtering and Smoothing", Cambridge
    %   University Press, 2013
    %
    % SEE ALSO
    %   GenericEKF, GenericGF, UnscentedTransform, GaussHermiteCubature
    %
    % VERSION
    %   2016-12-18
    % 
    % AUTHORS
    %   Roland Hostettler <roland.hostettler@aalto.fi>   
    
     %#ok<*PROPLC>
    
    %% Properties
    properties
        % The model of the type 'AGModel'
        model;
        
        % Mean and covariance
        m;
        P;
        
        % Predicted mean and covariance
        m_p;
        P_p;
                 
        % Stores the residual covariance (needed by IMMUKF)
        Pyy;
        v;
    end
    
    %% Private Properties
    properties (Access = private)
        % Sigma-point integration rule
        rule = UnscentedTransform();
    end

    %% Methods
    methods (Access = public)
        %% Constructor
        function self = AdditiveGF(model, rule)
            % Subclasses might not provide a model (e.g. the UKFIMM) as
            % they manage the model themselves. Hence, we have to check if
            % the model was passed or not.
            if nargin >= 1 && ~isempty(model)
                % Store the model
                self.model = model;
                self.initialize();
            end
            
            if nargin >= 2 && ~isempty(rule)
                self.rule = rule;
            end
        end
        
        %% Initialize the Filter
        function initialize(self)
            self.m = self.model.m0;
            self.P = self.model.P0;
        end
        
        %% Update the filter
        function update(self, y, t, u)
            self.timeUpdate(t, u);
            self.measurementUpdate(y, t, u);
        end
        
        %% Time update function
        function timeUpdate(self, t, u)
            % Augment the state and covariance
            m = self.m;
            P = self.P;
            Q = self.model.Q(m, t, u);
            Nx = size(m, 1);
            Nv = size(Q, 1);
            
            % Calculate & propagate the sigma points
            [X, wm, wc] = self.rule.calculateSigmaPoints(m, P);
            J = size(X, 2);
            X_p = zeros(Nx, J);
            for j = 1:J
                X_p(:, j) = self.model.f(X(1:Nx, j), zeros(Nv, 1), t, u);
            end
                
            % Predict the mean & covariance
            self.m_p = X_p*wm';
            P_p = zeros(Nx);
            for j = 1:J
                P_p = P_p + wc(j)*(X_p(:, j)-self.m_p)*(X_p(:, j)-self.m_p)';
            end
            P_p = P_p + Q;
            self.P_p = (P_p + P_p')/2;
        end
        
        %% Measurement update function
        function measurementUpdate(self, y, t, u)
            % Augment the predicted state and measurement covariance
            m_p = self.m_p;
            P_p = self.P_p;
            R = self.model.R(m_p, t, u);
            Nx = size(m_p, 1);
            Ny = size(y, 1);
            Nn = size(R, 1);

            % Propagate the sigma points through the measurement
            [X_p, wm, wc] = self.rule.calculateSigmaPoints(m_p, P_p);
            J = size(X_p, 2);
            Y_p = zeros(Ny, J);
            for j = 1:J
                Y_p(:, j) = self.model.g(X_p(1:Nx, j), zeros(Nn, 1), t, u);
            end

            % Predict the measurement and covariance
            y_p = Y_p*wm';
            Pyy = zeros(Ny);
            Pxy = zeros(Nx, Ny);
            for j = 1:J
                Pyy = Pyy + wc(j)*(Y_p(:, j) - y_p)*(Y_p(:, j) - y_p)';
                Pxy = Pxy ...
                    + wc(j)*(X_p(1:Nx, j) - self.m_p)*(Y_p(:, j) - y_p)';
            end
            Pyy = Pyy + R;
            Pyy = (Pyy+Pyy')/2;

            % Correction
            K = Pxy/Pyy;
            self.m = self.m_p + K*(y-y_p);
            Lyy = chol(Pyy, 'lower');
            P = self.P_p - (K*Lyy)*(K*Lyy)';
%             P = self.P_p - K*Pyy*K';

            % Ensure P is Hermitian
            self.P = (P + P')/2;
            self.Pyy = Pyy;
            self.v = y-y_p;
        end
    end
end
