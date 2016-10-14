classdef GenericCLGSSModel < CLGSSModel
    % Generic Mixed linear/non-linear Gaussian State Space Model
    % 
    % DESCRIPTION
    %   A generic mixed linear/non-linear Gaussian state space model of the
    %   form
    %
    %       xn[t] = fn(xn[t-1], u[t], t)
    %                   + An(xn[t-1], u[t], t) xl[t-1] + qn[t]
    %       xl[t] = fl(xn[t-1], u[t], t)
    %                   + Al(xn[t-1], u[t], t) xl[t-1] + ql[t]
    %       y[t] = h(xn[t], u[t], t) + C(xn[t], u[t], t) xl[t-1] + r[t]
    %
    %   where
    %
    %       * q[t] = [qn[t], ql[t]]^T ~ N(0, Q(xn[t-1], u[t], t),
    %       * r[t] ~ N(0, R(xn[t], u[t], t)),
    %       * x[0] = [xn[0], xl[0]]^T ~ N(m0, P0).
    %
    %   This class is essentially a convenience class for simple models
    %   where the different functions and/or matrices are easily defined in
    %   a script. More complex models should derive from the CLGSSModel
    %   base class directly and implement the corresponding function.
    %
    %   This model does not provide Jacobians.
    %
    % PROPERTIES
    %   None
    %
    % METHODS
    %   For a complete list of methods, see the parent class CLGSSModel.
    %
    %   Subclass-specific methods:
    %
    %   GenericCLGSSModel(fn, An, fl, Al, Q, h, C, R, m0, P0, in, il)
    %       Constructor that initializes the model. Requires all arguments.
    %
    %       fn, An, fl, Al, Q, h, C, R 
    %           Function handles of the form @(xn, u, t) corresponding to
    %           the respective function / matrix in the model above
    %
    %       m0, P0
    %           Initial state distribution mean and covariance
    %
    %       in, il
    %           Indices of the non-linear and linear states in the state
    %           vector such that
    %
    %               xn = x(in)
    %               xl = x(il)
    %
    %           where x is the overall state vector.
    %
    % SEE ALSO
    %   CLGSSModel, RBGF
    %
    % VERSION
    %   2016-10-14
    % 
    % AUTHORS
    %   Roland Hostettler <roland.hostettler@aalto.fi>   
    
    %% Properties
    properties (Access = private)
        % Initial state
        %m0;
        %P0;
        
        % Dynamics
        m_fn = @(xn, u, t) [];
        m_fl = @(xn, u, t) [];
        m_An = @(xn, u, t) [];
        m_Al = @(xn, u, t) [];
        m_Q = @(xn, u, t) [];
        
        % Observation
        m_h = @(xn, u, t) [];
        m_C = @(xn, u, t) [];
        m_R = @(xn, u, t) [];
    end

    %% Methods
    methods (Access = public)
        function self = GenericCLGSSModel(varargin)
            % Input order is: fn, An, fl, Al, Q, h, C, R, m0, P0, in, il
            if nargin < 12
                error('Model not fully specified.');
            end
            
            self.m_fn = varargin{1};
            self.m_An = varargin{2};
            self.m_fl = varargin{3};
            self.m_Al = varargin{4};
            self.m_Q = varargin{5};
            self.m_h = varargin{6};
            self.m_C = varargin{7};
            self.m_R = varargin{8};
            self.m0 = varargin{9};
            self.P0 = varargin{10};
            self.in = varargin{11};
            self.il = varargin{12};
        end
        
        function fn = fn(self, x, u, t)
            fn = self.m_fn(x, u, t);
        end
        
        function An = An(self, x, u, t)
            An = self.m_An(x, u, t);
        end
        
        function fl = fl(self, x, u, t)
            fl = self.m_fl(x, u, t);
        end
        
        function Al = Al(self, x, u, t)
            Al = self.m_Al(x, u, t);
        end

        function Q = Q(self, x, u, t)
            Q = self.m_Q(x, u, t);
        end
        
        function Qn = Qn(self)
            Q = self.m_Q(x, u, t);
            Qn = Q(self.in, self.in);
        end
        
        function Ql = Ql(self)
            Q = self.m_Q(x, u, t);
            Ql = Q(self.il, self.il);
        end
        
        function Qnl = Qnl(self)
            Q = self.m_Q(x, u, t);
            Qnl = Q(self.in, self.il);
        end
        
        function h = h(self, x, u, t)
            h = self.m_h(x, u, t);
        end
        
        function C = C(self, x, u, t)
            C = self.m_C(x, u, t);
        end
        
        function R = R(self, x, u, t)
            R = self.m_R(x, u, t);
        end
    end
end
