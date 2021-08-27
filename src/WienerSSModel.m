classdef WienerSSModel < GenericModel
    % Generic Wiener state-space model class
    % 
    % DESCRIPTION
    %   This is a generic class for Wiener-type state-space models, that is
    %   models with linear state dynamics and non-linear observations of
    %   the form
    %
    %       x[k] = A*x[k-1] + q[k]
    %       y[k] = g(x[k], r[k])
    %       x[0] ~ N(mu0, Sigma0)
    %
    %   The class is provides both a template for more specific 
    %   implementations but can also be used as a class on its own setting 
    %   the appropriate functions. Complex models should be children of
    %   this class since the generic class does not provide Jacobians, etc.
    %   However, this class it provides some convenience functions for easy
    %   subclassing (see below).
    %
    % PROPERTIES
    %   m0, P0
    %       The mean and covariance of the initial state
    %
    % METHODS
    %   WienerSSModel()
    %       Constructor to be used when subclassing (implicit).
    %
    %   WienerSSModel(F, Q, g, py_eval, mu0, Sigma0)
    %       Constructor to be used when using the generic class for a
    %       simple model where F and Q are time-invariant matrices,
    %       g(x, r, t, u) is the observation function, and 
    %       py_eval(y, x, t, u) the measurement likelihood. mu0 and Sigma0
    %       are the mean and covariance of the initial state.
    %
    %   F(t, u), Q(t, u)
    %       Convenience functions for the matrices A, Q, and R if these are
    %       time- and input-invariant. Subclasses can set the values for
    %       the corresponding matrix by accessing the protected m_A, m_Q,
    %       and m_R properties, respectively.
    %
    %   x = px0_rand(M)
    %       Random initial state generator, generates M random variables 
    %       such that x ~ N(mu0, Sigma0).
    %
    %   [x, Fx, Fq] = f(x, q, t, u)
    %       State transition function as defined above, valid for all
    %       Wiener-type systems. Returns the next state as well as the
    %       Jacobians w.r.t. x and q.
    %
    %   x = px_rand(x, q, t, u)
    %       Transition denisty random state generator, generates random
    %       variables according to x ~ p(x[k] | x[k-1]).
    %
    % SEE ALSO
    %   LGSSModel, WienerAPF, WienerAPS
    %
    % VERSION
    %   2016-10-18
    % 
    % AUTHORS
    %   Roland Hostettler <roland.hostettler@aalto.fi>
    
    % TODO
    %   * Include px_eval
    %   * Is the constructor really needed?
    
    %% Properties
    properties (Access = public)
        m0;
        P0;
    end
    
    %% 
    properties (Access = protected)
        m_F;
        m_Q;
        m_R;
    end

    %% 
    properties (Access = private)
        m_g;
        m_py_eval;
    end
    
    %% Public Methods
    methods (Access = public)
        %% Constructor
        function self = WienerSSModel(varargin)
            if nargin == 6
                self.m_F = varargin{1};
                self.m_Q = varargin{2};
                self.m_g = varargin{3};
                self.m_py_eval = varargin{4};
                self.m0 = varargin{5};
                self.P0 = varargin{6};
            elseif nargin ~= 0
                error('Model initialization failed, check the parameters.');
            end
        end
        
        %% Templated Functions
        % These are pre-implemented functions for time-invariant, static
        % values of A, Q, and R. Overwrite the functions as needed by the
        % model.
        function F = F(self, t, u)
            F = self.m_F;
        end
        
%         function B = B(self, t, u)
%             B = self.m_B;
%         end
                
        function Q = Q(self, x, t, u)
            Q = self.m_Q;
        end
        
        function yp = g(self, x, r, t, u)
            yp = self.m_g(x, r, t, u);
        end
        
        function py = py_eval(self, y, x, t, u)
            py = self.m_py_eval(y, x, t, u);
        end
        
        function R = R(self, x, t, u)
            R = self.m_R;
        end

        %% Generic Wiener State-Space System Functions
        % These are the same for all Wiener-type state-space systems and
        % ought not to be overwritten.
        function x = px0_rand(self, M)
            x = self.mu0*ones(1, M) ...
                + chol(self.Sigma0).'*randn(size(self.Sigma0, 1), M);
        end
        
        function [xp, Fx, Fq] = f(self, x, q, t, u)
%             xp = self.A(u, t)*x + self.B(u, t)*u + q;
            xp = self.F(t, u)*x + q;
            Fx = self.F(t, u);
            Fq = eye(size(q, 1));
        end
        
        function xp = px_rand(self, x, t, u)
            Q = self.Q(t, u);
            q = chol(Q).'*randn(size(Q, 1), size(x, 2));
            xp = self.f(x, q, t, u);
        end
        
        function px = px_eval(self, xp, x, t, u)
            Q = self.Q(t, u);
            px = mvnpdf( ...
                xp.', ...
                (self.f(x, zeros(size(Q, 1), size(x, 2)), t, u)).', ...
                self.Q ...
            ).';
        end
    end
end
