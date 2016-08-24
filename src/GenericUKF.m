classdef GenericUKF < handle
    % A generic unscented Kalman filter
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
    %   GenericEKF
    %
    % VERSION
    %   2016-03-02
    % 
    % AUTHORS
    %   Roland Hostettler <roland.hostettler@aalto.fi>   
    
    %% Properties
    properties
        % The model of the type 'GenericStateSpaceModel'
        model;
        
        % Mean and covariance
        m;
        P;
        
        % Predicted mean and covariance
        m_p;
        P_p;
        
        % Weight parameters
        alpha = 1e-3;
        beta = 2;
                 
        % Stores the residual covariance (needed by IMMUKF)
        Pyy;
        v;
    end

    %% Methods
    methods (Access = public)
        %% Constructor
        function self = GenericUKF(model)
            % Subclasses might not provide a model (e.g. the UKFIMM) as
            % they manage the model themselves. Hence, we have to check if
            % the model was passed or not.
            if nargin == 1
                % Store the model
                self.model = model;

                % Initialize the filter
                self.m = model.m0;
                self.P = model.P0;
            end
        end
        
        %% Update the filter
        function update(self, u, y, t)
            self.timeUpdate(u, t);
            self.measurementUpdate(u, y, t);
            self.model.t = t;
        end
        
        %% Time update function
        function timeUpdate(self, u, t)
            % Augment the state and covariance
            Q = self.model.Q;
            Nx = size(self.m, 1);
            Nv = size(Q, 1);
            x_a = [self.m; zeros(Nv, 1)];
            P_a = [ ...
                       self.P, zeros(Nx, Nv); ...
                zeros(Nv, Nx),             Q; ...
            ];            
            L = Nx+Nv;
            lambda = (self.alpha^2-1)*L;
            [Wm, Wc] = self.weights(L);
            
            % Calculate the sigma points
            dx = sqrt(L+lambda)*chol(P_a, 'lower');
            X_a = [x_a, x_a*ones(1, L) + dx, x_a*ones(1, L) - dx];
            
            % Propagate all the sigma points
            X_p = zeros(Nx, 2*L+1);
            for n = 1:2*L+1
                X_p(:, n) = self.model.f( ...
                    X_a(1:Nx, n), ...
                    u, ...
                    X_a(Nx+1:Nx+Nv, n), ...
                    t ...
                );
            end
                
            % Predict the mean & covariance
            self.m_p = X_p*Wm';
            P_p = zeros(Nx);
            for n = 1:2*L+1
                P_p = P_p + Wc(n)*(X_p(:, n)-self.m_p)*(X_p(:, n)-self.m_p)';
            end
            
            % Ensure P_p is Hermitian
%             self.P_p = P_p;
            self.P_p = (P_p + P_p')/2;
        end
        
        %% Measurement update function
        function measurementUpdate(self, u, y, t)
                % Augment the predicted state and measurement covariance
                R = self.model.R;
                Nx = size(self.m, 1);
                Ny = size(y, 1);
                Nn = size(R, 1);
                x_p = [self.m_p; zeros(Nn, 1)];
                P_a = [ ...
                         self.P_p, zeros(Nx, Nn); ...
                    zeros(Nn, Nx),             R; ...
                ];
                L = Nx + Nn;
                lambda = (self.alpha^2-1)*L;
                [Wm, Wc] = self.weights(L);

                % Calculate the sigma points
                dx = sqrt(L+lambda)*chol(P_a, 'lower');
                X_p = [x_p, x_p*ones(1, L) + dx, x_p*ones(1, L) - dx];

                % Propagate the sigma points through the measurement
                Y_p = zeros(Ny, 2*L+1);
                for n = 1:2*L+1
                    Y_p(:, n) = self.model.g(...
                        X_p(1:Nx, n), ...
                        u, ...
                        X_p(Nx+1:Nx+Nn, n), ...
                        t ...
                    );
                end

                % Predict the measurement and covariance
                y_p = Y_p*Wm';
                Pyy = zeros(Ny);
                Pxy = zeros(Nx, Ny);
                for n = 1:2*L+1
                    Pyy = Pyy + Wc(n)*(Y_p(:, n) - y_p)*(Y_p(:, n) - y_p)';
                    Pxy = Pxy ...
                        + Wc(n)*(X_p(1:Nx, n) - self.m_p)*(Y_p(:, n) - y_p)';
                end

                % Ensure Pyy is Hermitian
                Pyy = (Pyy+Pyy')/2;

                % Correction
                K = Pxy/Pyy;
                self.m = self.m_p + K*(y-y_p);
                Lyy = chol(Pyy, 'lower');
                P = self.P_p - (K*Lyy)*(K*Lyy)';
%                 P = self.P_p - K*Pyy*K';

                % Ensure P is Hermitian
%                 self.P = P;
                self.P = (P + P')/2;
                
                % Store additional items
                self.Pyy = Pyy;
                self.v = y-y_p;
        end
       
        %% Calculates the mean- and covariance weights
        function [Wm, Wc] = weights(self, L)
            lambda = (self.alpha^2-1)*L;
            Wm = 1/(2*(L+lambda))*ones(1, 2*L+1);
            Wm(1) = lambda/(L+lambda);
            Wc = Wm;
            Wc(1) = lambda/(L+lambda) + (1-self.alpha^2+self.beta);
        end
    end
end
