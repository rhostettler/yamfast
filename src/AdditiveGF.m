classdef AdditiveGF < GaussianFilter
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
    %   Does not provide properties of its own, see the GaussianFilter's
    %   documentation for inherited properties.
    %
    % METHODS
    %   Implements all the function as required (and described) by the
    %   GaussianFilter-class. See its documentation for details.
    %
    %   AdditiveGF(model, rule)
    %       Constructor to initialize the filter.
    %
    %       model (optional)
    %           State-space model for this filter to operate on. While
    %           specifying a model is optional, omitting the model is only
    %           meaningful in cases where this filter is used together with
    %           an IMM (see IMMFilter). In all other cases, the model must
    %           be specified.
    %
    %       rule (optional)
    %           Sigma-point method to use. If not specified, the unscented
    %           transform is used with its default parameters (see
    %           UnscentedTransform for details), effectively giving an UKF.
    %   
    % REFERENCES
    %   [1] S. S?rkk?, "Bayesian Filtering and Smoothing", Cambridge
    %       University Press, 2013
    %
    % SEE ALSO
    %   GaussianFilter, KalmanFilter, ExtendedKF, GenericGF
    %
    % VERSION
    %   2016-12-22
    % 
    % AUTHORS
    %   Roland Hostettler <roland.hostettler@aalto.fi>   
    
     %#ok<*PROPLC>
    
    %% Properties
    properties
        % No extra properties
    end
    
    %% Private Properties
    properties (Access = private)
        % Sigma-point integration rule, UT by default
        rule = UnscentedTransform();
    end

    %% Methods
    methods (Access = public)
        %% Constructor
        function self = AdditiveGF(model, rule)
            if nargin >= 1 && ~isempty(model)
                self.model = model;
                self.initialize();
            end
            if nargin >= 2 && ~isempty(rule)
                self.rule = rule;
            end
        end
                
        %% Time Update
        function [m_p, P_p]  = timeUpdate(self, t, u)
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
            m_p = X_p*wm';
            P_p = zeros(Nx);
            for j = 1:J
                P_p = P_p + wc(j)*(X_p(:, j)-self.m_p)*(X_p(:, j)-self.m_p)';
            end
            P_p = P_p + Q;
            P_p = (P_p + P_p')/2;
            
            self.m_p = m_p;
            self.P_p = P_p;
        end
        
        %% Measurement Update
        function [m, P] = measurementUpdate(self, y, t, u)
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
            m = self.m_p + K*(y-y_p);
            Lyy = chol(Pyy, 'lower');
            P = self.P_p - (K*Lyy)*(K*Lyy)';
%             P = self.P_p - K*Pyy*K';
            P = (P + P')/2;

            self.m = m;
            self.P = P;
            self.Pyy = Pyy;
            self.v = y-y_p;
        end
    end
end
