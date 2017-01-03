classdef AWGNModel < GenericModel
    % Additive white Gaussian noise model
    % 
    % DESCRIPTION
    %   Model for non-linear state-space systems with additive, white
    %   Gaussian process- and measurement noise of the form
    %
    %       x[n] = f(x[n-1], t[n], u[n]) + q[n]
    %       y[n] = g(x[n], t[n], u[n]) + r[n]
    %
    %   where q[n] ~ N(0, Q), r[n] ~ N(0, R), and x[0] ~ N(m0, P0).
    %
    %   The class can be used in two ways:
    %       
    %       1. As its own class for simple models where the state
    %          transition and measurement functions can be described by
    %          function handles and the process and measurement noise
    %          matrices are constant.
    %
    %       2. As the parent class for more complex, custom models. In this
    %          case, the child class should re-implement the functions f,
    %          Q, g, and R.
    %
    % PROPERTIES
    %   m0, P0
    %       Mean and covaraince of the initial state.
    %
    % METHODS
    %   AWGNModel(f, Fx, Q, g, Gx, R)
    %       Constructor to be used when the model is used in its simple
    %       form. Then, f, Q, g, and R have to be specified and must be of
    %       the following form:
    %
    %       f = @(x, t, u)
    %           Function handle for the state transition function.
    %
    %       Fx = @(x, t, u)
    %           Function handle for the Jacobian of the state transition
    %           function. May be left empty (but required if EKF-type
    %           filters are used).
    %
    %       Q   Process noise covariance matrix.
    %       
    %       g = @(x, t, u)
    %           Function handle for the measurement function.
    %
    %       Gx = @(x, t, u)
    %           Function handle for the Jacobian of the measurement
    %           function. May be left empty (but required if EKF-type
    %           filters are used).
    %
    %       R   Measurement noise covariance matrix.
    %
    % SEE ALSO
    %   GenericModel
    %
    % VERSION
    %   2016-12-22
    % 
    % AUTHORS
    %   Roland Hostettler <roland.hostettler@aalto.fi>   
    
    %% Public Properties
    properties (Access = public)
        % Mean and covariance of the initial state
        m0;
        P0;
    end

    %% Protected Properties
    properties (Access = private)
        m_f;
        m_Fx;
        m_Q;
        m_g;
        m_Gx;
        m_R;
    end
            
    %% Common Methods
    methods (Access = public)
        %% Constructor
        % Only used when the 'simple' form of the model is used.
        function self = AWGNModel(f, Fx, Q, g, Gx, R)
            switch nargin
                case 0
                    % nop
                    
                case 6
                    self.m_f = f;
                    self.m_Fx = Fx;
                    self.m_Q = Q;
                    self.m_g = g;
                    self.m_Gx = Gx;
                    self.m_R = R;
                    
                otherwise
                    error('Model not fully specified, see model documentation.');
            end
        end
        
        %% Draw Samples from the Initial Distribution
        function x0 = px0_rand(self, M)
            x0 = self.m0*ones(1, M) ...
                + chol(self.P0).'*randn(size(self.P0, 1), M);
        end
        
        %% Default State Transition Function
        function [xp, Fx, Fq] = f(self, x, q, t, u)
            if ~isempty(self.m_Fx)
                Fx = self.m_Fx(x, t, u);
            else
                Fx = [];
            end
            Fq = eye(size(q, 1));
            xp = self.m_f(x, t, u) + q;
        end
        
        %% Default Process Noise Covariance
        function Q = Q(self, x, t, u)
            Q = self.m_Q;
        end
        
        %% Draw Samples from the State Transition Density
        function xp = px_rand(self, x, t, u)
            Q = self.Q(x, t, u);
            q = chol(Q).'*randn(size(Q, 1), size(x, 2));
            xp = self.f(x, q, t, u);
        end
        
        %% Evaluate the State Transition Density
        function px = px_eval(self, xp, x, t, u)
            Q = self.Q(x, t, u);
            f = self.f(x, t, zeros(size(Q, 1), size(x, 2)), u);
            px = mvnpdf(xp.', f.', Q).';
        end
        
        %% Default Measurement Function
        function [y, Gx, Gr] = g(self, x, r, t, u)
            if ~isempty(self.m_Gx)
                Gx = self.m_Gx(x, t, u);
            else
                Gx = [];
            end
            Gr = eye(size(r, 1));
            y = self.g(x, t, u) + r;
        end
        
        %% Default Measurement Noise Covariance
        function R = R(self, x, t, u)
            R = self.m_R;
        end
                
        %% Evaluate the Likelihood
        function py = py_eval(self, y, x, t, u)
            R = self.R(x, t, u);
            py = mvnpdf(y.', self.g(x, t, u).', R).';
        end
    end
end
