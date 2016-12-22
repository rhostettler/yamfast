classdef GenericGF < GaussianFilter
    % Generic Gaussian assumed density sigma-point filter
    % 
    % DESCRIPTION
    %   Gaussian assumed density filter for non-linear state space models
    %   with non-additive process- and measurement noise of the form
    %
    %       x[n] = f(x[n-1], q[n], t[n], u[n])
    %       y[n] = g(x[n], r[n], t[n], u[n])
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
    %   GenericGF(model, rule)
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
    % SEE ALSO
    %   GaussianFilter, AdditiveGF
    %
    % VERSION
    %   2016-12-22
    % 
    % AUTHORS
    %   Roland Hostettler <roland.hostettler@aalto.fi>
        
    %% Properties
    properties
        % No properties of its own
    end
    
    %% Private Properties
    properties (Access = private)
        % Sigma-point integration rule, UT by default
        rule = UnscentedTransform();
    end
    
    %% Methods
    methods (Access = public)
        %% Constructor
        function self = GenericGF(model, rule)
            if nargin >= 1 && ~isempty(model)
                self.model = model;
                self.initialize();
            end
            if nargin >= 2 && ~isempty(rule)
                self.rule = rule;
            end
        end
        
        %% Time Update
        function [m_p, P_p] = timeUpdate(self, t, u)
            % Augment the state and covariance
            Q = self.model.Q(self.m, t, u);
            Nx = size(self.m, 1);
            Nq = size(Q, 1);
            x_a = [self.m; zeros(Nq, 1)];
            P_a = [
                       self.P, zeros(Nx, Nq);
                zeros(Nq, Nx),             Q;
            ];
        
            % Calculate & propagate the sigma points
            [X_a, wm, wc] = self.rule.calculateSigmaPoints(x_a, P_a);
            L = size(X_a, 2);
            X_p = zeros(Nx, L);
            for l = 1:L
                X_p(:, l) = self.model.f(X_a(1:Nx, l), X_a(Nx+1:Nx+Nq, l), t, u);
            end
                
            % Predict the mean & covariance
            m_p = X_p*wm';
            P_p = zeros(Nx);
            for l = 1:L
                P_p = P_p + wc(l)*(X_p(:, l)-self.m_p)*(X_p(:, l)-self.m_p)';
            end
            P_p = (P_p + P_p')/2;
%             self.P_p = self.stabilize(P_p);

            self.m_p = m_p;
            self.P_p = P_p;
        end
        
        %% Measurement update function
        function [m, P] = measurementUpdate(self, y, t, u)
                R = self.model.R(self.m_p, t, u);
                Nx = size(self.m, 1);
                Ny = size(y, 1);
                Nr = size(R, 1);
                x_a = [self.m_p; zeros(Nr, 1)];
                P_a = [
                         self.P_p, zeros(Nx, Nr);
                    zeros(Nr, Nx),             R;
                ];
            
                % Calculate and propagate sigma points
                [X_p, Wm, Wc] = self.rule.calculateSigmaPoints(x_a, P_a);
                L = size(X_p, 2);
                Y_p = zeros(Ny, 2*L+1);
                for l = 1:L
                    Y_p(:, l) = self.model.g(X_p(1:Nx, l), X_p(Nx+1:Nx+Nr, l), t, u);
                end

                % Predict the measurement and covariance
                y_p = Y_p*Wm';
                Pyy = zeros(Ny);
                Pxy = zeros(Nx, Ny);
                for l = 1:2*L+1
                    Pyy = Pyy + Wc(l)*(Y_p(:, l) - y_p)*(Y_p(:, l) - y_p)';
                    Pxy = Pxy ...
                        + Wc(l)*(X_p(1:Nx, l) - self.m_p)*(Y_p(:, l) - y_p)';
                end

                % Ensure Pyy is Hermitian
                Pyy = (Pyy+Pyy')/2;

                % Correction
                K = Pxy/Pyy;
                m = self.m_p + K*(y-y_p);
                P = self.P_p - K*Pyy*K';
                P = (P + P')/2;
                %P = self.stabilize(P);

                self.m = m;
                self.P = P;
                self.Pyy = Pyy;
                self.v = y-y_p;
        end
    end

    %% Private Methods
    methods (Access = private)
        function P = stabilize(self, P)
            epsilon = 1e-15;
            P = (P+P')/2;
            [U, V] = eig(P);
            V = diag(V);
            if min(V) < epsilon
                V(V < epsilon) = epsilon;
%                 P = P+2*abs(lmin)*eye(size(P, 1));
                P = U*diag(V)*U';
                warning('Covariance stabilized.');
            end
        end
    end
end
