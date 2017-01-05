classdef LGSSModel < AWGNModel
    % Linear, Gaussian state-space model
    % 
    % DESCRIPTION
    %   Linear, Gaussian state-space model of the form
    %
    %       x[n] = F(t[n], u[n]) x[n-1] + q[n]
    %       y[n] = G(t[n], u[n]) x[k] + r[n]
    %
    %   where
    %
    %       * q[n] ~ N(0, Q)
    %       * r[n] ~ N(0, R)
    %       * x[0] ~ N(m0, P0)
    %
    % PROPERTIES
    %   Inherits the properties from AWGNModel.
    %
    % METHODS
    %   LGSSModel(F, G, Q, R, m0, P0)
    %       Initialization of the model where:
    %
    %           F   State transition matrix;
    %           G   Observation matrix;
    %           Q   Process noise covariance matrix;
    %           R   Measurement noise covariance matrix;
    %           m0  Mean of the initial state;
    %           P0  Covariance of the initial state.
    %
    %       All parameters can either be function handles of the form 
    %       @(t, u) or plain matrices.
    %
    % SEE ALSO
    %   AWGNModel, GenericModel
    %
    % VERSION
    %   2017-01-02
    % 
    % AUTHORS
    %   Roland Hostettler <roland.hostettler@aalto.fi>   
        
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
        function self = LGSSModel(F, G, Q, R, m0, P0)
            self.m_F = F;
            self.m_G = G;
            self.m_Q = Q;
            self.m_R = R;
            self.m0 = m0;
            self.P0 = P0;
        end
                
        %% Process Dynamics Matrix
        function F = F(self, t, u)
            F = self.get(self.m_F, t, u);
        end
        
        %% Process Noise Covariance
        function Q = Q(self, x, t, u)
            Q = self.get(self.m_Q, t, u);
        end
        
        %% State Transition Function
        function [xp, Fx, Fq] = f(self, x, q, t, u)
            Fx = self.F(t, u);
            xp = Fx*x + q;
            Fq = eye(size(q, 1));
        end
                
        %% Observation Matrix
        function G = G(self, t, u)
            G = self.get(self.m_G, t, u);
        end
        
        %% Measurement Noise Covariance
        function R = R(self, x, t, u)
            R = self.get(self.m_R, t, u);
        end
        
        %% Measurement Function
        function [y, Gx, Gr] = g(self, x, r, t, u)
            Gx = self.G(t, u);
            y = self.m_G*x + r;
            Gr = eye(size(r, 1));
        end
    end
    
    %% Private Helper Methods
    methods (Access = private)
        function Y = get(self, X, t, u)
            Y = X;
            if isa(X, 'function_handle')
                Y = X(t, u);
            end
        end
    end
end
