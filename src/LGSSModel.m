classdef LGSSModel < GenericModel
    % Linear, Gaussian state-space model
    % 
    % DESCRIPTION
    %   Linear, Gaussian state-space model of the form
    %
    %       x(k) = F(t(k)) x(k-1) + B(t(k)) u(k-1) q(k-1)
    %       y(k) = G(t(k)) x(k) + D(t(k)) u(k) + r(k)
    %
    %   where
    %
    %       * q(k) ~ N(0, Q)
    %       * r(k) ~ N(0, R)
    %       * x(0) ~ N(mu0, Sigma0)
    %
    % PROPERTIES
    %   mu0
    %       Mean of the initial state.
    %
    %   Sigma0
    %       Covariance of the initial state.
    %
    % METHODS
    %   LGSSModel(F, G, Q, R, mu0, Sigma0)
    %       Initialization of the model.
    %
    % SEE ALSO
    %   GenericModel
    %
    % VERSION
    %   2016-10-19
    % 
    % AUTHORS
    %   Roland Hostettler <roland.hostettler@aalto.fi>   
    
    %% Properties
    properties (Access = public)
        % Initial state and covariance
        mu0;
        Sigma0;
    end
    
    %% Private Properties
    properties (Access = private)
        m_F;
        m_Q;
        m_G;
        m_R;
    end

    %% Methods
    methods (Access = public)
        %% Constructor
        function self = LGSSModel(F, G, Q, R, mu0, Sigma0)
            self.m_F = F;
            self.m_G = G;
            self.m_Q = Q;
            self.m_R = R;
            self.mu0 = mu0;
            self.Sigma0 = Sigma0;
        end
        
        %% Initial Random State Generator
        function x0 = px0_rand(self, M)
            x0 = self.mu0*ones(1, M) ...
                + chol(self.Sigma0).'*randn(size(self.Sigma0, 1), M);
        end
        
        %% Process Dynamics Matrix
        function F = F(self, t, u)
            F = self.m_F;
        end
        
        %% Process Noise Covariance
        function Q = Q(self, t, u)
            Q = self.m_Q;
        end
        
        %% State Transition Function
        function [xp, Fx, Fq] = f(self, x, q, t, u)
%             xp = self.m_F*x + self.m_B*u + q;
            xp = self.m_F*x + q;
            Fx = self.m_F;
            Fq = eye(size(q, 1));
        end
        
        %% State Transition Random State Generator
        function xp = px_rand(self, x, t, u)
            q = chol(self.m_Q).'*randn(size(self.m_Q, 1), size(x, 2));
            xp = self.f(x, q, t, u);
        end
        
        %% Observation Matrix
        function G = G(self, t, u)
            G = self.m_G;
        end
        
        %% Measurement Noise Covariance
        function R = R(self, t, u)
            R = self.m_R;
        end
        
        %% Measurement Function
        function [y_p, Gx, Gr] = g(self, x, r, t, u)
%             M = size(x, 2);
%             y_p = self.m_G*x + self.m_D*u*ones(1, M) + r;
            y_p = self.m_G*x + r;
            Gx = self.m_G;
            Gr = eye(size(r, 1));
        end
        
        %% Likelihood
        function py = py_eval(self, y, x, t, u)
            py = mvnpdf(y.', (self.m_G*x).', self.m_R).';
        end
    end
end
