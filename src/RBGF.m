classdef RBGF < handle
    % Rao-Blackwellized Gaussian Filter
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
    %   2016-08-25
    % 
    % AUTHORS
    %   Roland Hostettler <roland.hostettler@aalto.fi>   
    
    % TODO
    %   * Generalize to a general Gaussian filter
    %   * Implement posterior linearization update
    %   * Implement alias properties for mn, ml, Pn, etc.
    %   * Implement for hierarchical model
    
    % Properties
    properties
        model = [];
        
        % State
        m;
        P;
        
        % Predicted state
        m_p;
        P_p;
        
        % UKF parameters
        alpha = 1e-3;
        beta = 2;
        kappa = 0;
    end

    % Methods
    methods
        function self = RBGF(model)
            if nargin == 1
                self.model = model;
                
                self.m = model.m0;
                self.P = model.P0;
            end
        end
        
        function update(self, u, y, t)
            self.timeUpdate(u, t);
            self.measurementUpdate(u, y, t);
            self.model.t = t;
        end
        
        %% 
        function timeUpdate(self, u, t)
            model = self.model;
            in = model.in;
            il = model.il;
            
            xn = self.m(in);
            Nxn = size(xn, 1);
            xl = self.m(il);
            Nxl = size(xl, 1);
            Nx = Nxn+Nxl;
            Nq = size(model.Q(xn, u, t), 1);
            
            P = self.P;
            Pn = self.P(in, in);
            Pl = self.P(il, il);
            Pnl = self.P(in, il);
            
            % Calculate the sigma-points
            [Xn, wm, wc] = self.calculateSigmaPoints(xn, Pn);
            M = size(Xn, 2);
            
            % Whiten
            L = Pnl'/Pn;
            Pl_tilde = Pl - L*Pn*L';
            Xl_tilde = xl*ones(1, M) + L*(Xn - xn*ones(1, M));
            
            % Propagate the sigma points
            X = zeros(Nx, M);
            X_p = zeros(Nx, M);
            for m = 1:M
                X(in, m) = Xn(:, m);
                X(il, m) = Xl_tilde(:, m);
                X_p(:, m) = model.f(X(:, m), u, zeros(Nq, 1), t);
            end
            
            % Calculate the mean and covariances
            x_p = X_p*wm(:);
            P_p = zeros(Nx, Nx);
            A = zeros(Nx, Nxl);
            for m = 1:M
                A(in, :) = model.An(Xn(:, m), u, t);
                A(il, :) = model.Al(Xn(:, m), u, t);
                Q = model.Q(Xn(:, m), u, t);
                P_p = P_p + wc(m)*( ...
                    (X_p(:, m) - x_p)*(X_p(:, m) - x_p)' + A*Pl_tilde*A' + Q ...
                );
            end            
            
            % Store
            self.m_p = x_p;
            self.P_p = P_p;
        end
        
        %% 
        function measurementUpdate(self, u, y, t)
            model = self.model;
            in = model.in;
            il = model.il;
            x_p = self.m_p;
            xn_p = x_p(in);
            Nxn = size(xn_p, 1);
            xl_p = x_p(il);
            Nxl = size(xl_p, 1);
            Nx = Nxl+Nxn;
            P_p = self.P_p;
            Pn_p = P_p(in, in);
            Pnl_p = P_p(in, il);
            Pl_p = P_p(il, il);
            Ny = size(y, 1);
            
            % Calculate the sigma-points
            [Xn_p, wm, wc] = calculateSigmaPoints(self, xn_p, Pn_p);
            M = size(Xn_p, 2);
            
            % Whiten
            L = Pnl_p'/Pn_p;
            Xl_tilde_p = xl_p*ones(1, M) + L*(Xn_p - xn_p*ones(1, M));
            Pl_tilde_p = Pl_p - L*Pn_p*L';
            
            % Propagate the sigma-points
            X = zeros(Nx, M);
            Y_p = zeros(Ny, M);
            for m = 1:M
                X(in, m) = Xn_p(:, m);
                X(il, m) = Xl_tilde_p(:, m);
                Y_p(:, m) = model.g(X(:, m), u, zeros(size(model.R)), t);
            end
            
            % Calculate the predicted output and covariances
            y_p = Y_p*wm(:);
            Dn = zeros(Nxn, Ny);
            Dl = zeros(Nxl, Ny);
            S = zeros(Ny, Ny);
            for m = 1:M
                Dn = Dn + wc(m)*((Xn_p(:, m) - xn_p)*(Y_p(:, m) - y_p)');
                Dl = Dl + wc(m)*model.C(Xn_p(:, m), u, t)';
                S = S + wc(m)*( ...
                    (Y_p(:, m) - y_p)*(Y_p(:, m) - y_p)' ...
                    + model.C(Xn_p(:, m), u, t)*Pl_tilde_p*model.C(Xn_p(:, m), u, t)' ...
                    + model.R(Xn_p(:, m), u, t) ...
                );
            end
            Dl = Pl_tilde_p*Dl + L*Dn;
            
            % Update mean and covariance
            D = zeros(Nx, Ny);
            D(in, :) = Dn;
            D(il, :) = Dl;
            K = D/S;
            x = x_p + K*(y-y_p);
            P = P_p - K*S*K';
            
            % Store
            self.m = x;
            self.P = P;
        end
        
        %% Cubature Sigma Points
%         function [X, wm, wc] = calculateSigmaPoints(self, mu, Sigma)
%             N = size(mu, 1);
%             M = 2*N;
%             wm = 1/M*ones(1, M);
%             wc = 1/M*ones(1, M);
%             L = chol(Sigma, 'lower');
%             X = mu*ones(1, M) + sqrt(N)*[L -L];
%         end
        
        %% Unscented Sigma Points
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
