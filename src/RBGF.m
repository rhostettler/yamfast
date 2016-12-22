classdef RBGF < handle
    % Rao-Blackwellized Gaussian filter
    % 
    % DESCRIPTION
    %   Rao-Blackwellized Gaussian filter for conditionally linear Gaussian
    %   state space systems.
    %
    % PROPERTIES
    % 
    %
    % METHODS
    %   update(y, t, u)
    %
    % REFERENCES
    %   [1] R. Hostettler and S. Särkkä, "Rao-Blackwellized Gaussian
    %       Filtering and Smoothing", 2017.
    %
    % SEE ALSO
    %
    %
    % VERSION
    %   2016-12-18
    % 
    % AUTHORS
    %   Roland Hostettler <roland.hostettler@aalto.fi>   
    
    % TODO
    %   * Generalize to a general Gaussian filter
    %   * Implement alias properties for mn, ml, Pn, etc.
    %   * Implement for hierarchical model
    
    %#ok<*PROPLC>
    
    %% Properties
    properties (Access = public)
        model = [];
        
        % State estimate and covariance
        m;
        P;
        
        % Predicted state and covaraince
        m_p;
        P_p;
        
        % Number of measurement update iterations (= posterior
        % linearization approximation if not equal to 1)
        Ni = 1;
    end
    
    properties (Access = private)
        % The sigma-point rule
        rule = UnscentedTransform();
    end
    
    %% Property Setter/Getter Methods
    methods
        function model.set(self, model)
            self.model = model;
            self.initialize();
        end
    end
    
    %% Public Methods
    methods (Access = public)
        %% Constructor
        function self = RBGF(model, rule)
            if nargin >= 1 && ~isempty(model)
                self.model = model;
                self.initialize();
            end
            if nargin >= 2 && ~isempty(rule)
                self.rule = rule;
            end
        end
        
        %% Initialization
        function initialize(self)
            self.m = self.model.m0;
            self.P = self.model.P0;
        end
        
        %% Filter Update
        function update(self, y, t, u)
            self.timeUpdate(t, u);
            self.measurementUpdate(y, t, u);
        end
        
        %% Time Update (Prediction)
        function timeUpdate(self, t, u)
            model = self.model;
            in = model.in;
            il = model.il;
            mn = self.m(in);
            Nxn = size(mn, 1);
            ml = self.m(il);
            Nxl = size(ml, 1);
            Nx = Nxn+Nxl;
            Nq = size(model.Q(mn, t, u), 1);
            Pn = self.P(in, in);
            Pl = self.P(il, il);
            Pnl = self.P(in, il);
            
            % Calculate the sigma-points
            [Xn, wm, wc] = self.rule.calculateSigmaPoints(mn, Pn);
            J = size(Xn, 2);
            
            % Whiten
            L = Pnl'/Pn;
            Pl_tilde = Pl - L*Pn*L';
            Xl_tilde = ml*ones(1, J) + L*(Xn - mn*ones(1, J));
            
            % Propagate the sigma points
            X = zeros(Nx, J);
            X_p = zeros(Nx, J);
            for j = 1:J
                X(in, j) = Xn(:, j);
                X(il, j) = Xl_tilde(:, j);
                X_p(:, j) = model.f(X(:, j), zeros(Nq, 1), t, u);
            end
            
            % Calculate the mean and covariances
            m_p = X_p*wm(:);
            P_p = zeros(Nx, Nx);
            A = zeros(Nx, Nxl);
            Q = zeros(Nx, Nx);
            for j = 1:J
                A(in, :) = model.An(Xn(:, j), t, u);
                A(il, :) = model.Al(Xn(:, j), t, u);
                Qn = model.Qn(Xn(:, j), t, u);
                Ql = model.Ql(Xn(:, j), t, u);
                Qnl = model.Qnl(Xn(:, j), t, u);
                Q(in, in) = Qn;
                Q(in, il) = Qnl;
                Q(il, in) = Qnl';
                Q(il, il) = Ql;
                P_p = P_p + wc(j)*( ...
                    (X_p(:, j) - m_p)*(X_p(:, j) - m_p)' + A*Pl_tilde*A' + Q ...
                );
            end            
            
            % Store
            self.m_p = m_p;
            self.P_p = (P_p+P_p')/2;
        end
        
        %% Measurement Update
        function measurementUpdate(self, y, t, u)
            m_i = self.m_p;
            P_i = self.P_p;
            for i = 1:self.Ni
                [m_i, P_i] = self.measurementUpdateIteration(y, t, u, m_i, P_i);
            end
            self.m = m_i;
            self.P = P_i;
        end
    end
    
    %% Internal Methods
    methods (Access = protected)
        %% Single Measurement Update Iteration
        function [m, P] = measurementUpdateIteration(self, y, t, u, m_p, P_p)
            model = self.model;
            in = model.in;
            il = model.il;
            mn_p = m_p(in);
            Nxn = size(mn_p, 1);
            ml_p = m_p(il);
            Nxl = size(ml_p, 1);
            Nx = Nxl+Nxn;
            Pn_p = P_p(in, in);
            Pnl_p = P_p(in, il);
            Pl_p = P_p(il, il);
            Ny = size(y, 1);
            Nr = size(model.R(m_p, t, u), 1);
            
            % Calculate the sigma-points
            [Xn_p, wm, wc] = self.rule.calculateSigmaPoints(mn_p, Pn_p);
            J = size(Xn_p, 2);
            
            % Whiten
            L = Pnl_p'/Pn_p;
            Xl_tilde_p = ml_p*ones(1, J) + L*(Xn_p - mn_p*ones(1, J));
            Pl_tilde_p = Pl_p - L*Pn_p*L';
            
            % Propagate the sigma-points
            X = zeros(Nx, J);
            Y_p = zeros(Ny, J);
            for j = 1:J
                X(in, j) = Xn_p(:, j);
                X(il, j) = Xl_tilde_p(:, j);
                Y_p(:, j) = model.g(X(:, j), zeros(Nr, 1), t, u);
            end
            
            % Calculate the predicted output and covariances
            y_p = Y_p*wm(:);
            Dn = zeros(Nxn, Ny);
            Dl = zeros(Nxl, Ny);
            S = zeros(Ny, Ny);
            for j = 1:J
                Dn = Dn + wc(j)*((Xn_p(:, j) - mn_p)*(Y_p(:, j) - y_p)');
                Dl = Dl + wc(j)*model.C(Xn_p(:, j), t, u)';
                S = S + wc(j)*( ...
                    (Y_p(:, j) - y_p)*(Y_p(:, j) - y_p)' ...
                    + model.C(Xn_p(:, j), t, u)*Pl_tilde_p*model.C(Xn_p(:, j), t, u)' ...
                    + model.R(Xn_p(:, j), t, u) ...
                );
            end
            Dl = Pl_tilde_p*Dl + L*Dn;
            
            % Update mean and covariance
            D = zeros(Nx, Ny);
            D(in, :) = Dn;
            D(il, :) = Dl;
            K = D/S;
            m = m_p + K*(y-y_p);
            P = P_p - K*S*K';
            P = (P+P')/2;
        end
    end
end
