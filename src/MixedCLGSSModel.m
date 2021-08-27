classdef MixedCLGSSModel < AWGNModel
    % Mixed linear/non-linear Gaussian state-space model
    % 
    % DESCRIPTION
    %   A mixed linear/non-linear Gaussian state-space model of the form
    %
    %       xn[n] = fn(xn[n-1], u[n], t[n])
    %                   + An(xn[n-1], u[n], t[n]) xl[n-1] + qn[n]
    %       xl[n] = fl(xn[n-1], u[n], t[n])
    %                   + Al(xn[n-1], u[n], t[n]) xl[n-1] + ql[n]
    %       y[n] = h(xn[n], u[n], t[n]) + C(xn[n], u[n], t[n]) xl[n-1] 
    %                   + r[n]
    %   where
    %
    %       * q[n] = [qn[n], ql[n]]^T ~ N(0, Q(xn[n-1], u[n], t[n]),
    %       * r[n] ~ N(0, R(xn[n], u[n], t[n])),
    %       * x[0] = [xn[0], xl[0]]^T ~ N(m0, P0).
    %
    %   This class can be used in two ways:
    %
    %       1. As a convenience class for simple models where the different
    %          functions and/or matrices are easily defined as function 
    %          handles in a script.
    %
    %       2. As the parent class for more complex models that directly
    %          implement the corresponding functions. In this form, all the
    %          functions must be implemented by the subclass.
    %
    % PROPERTIES
    %   This class inhertis the properties defined by the AWGNModel and
    %   defines the following additional properties.
    %
    %   in (r/w, default: [])
    %       Vector of indices to the non-linear states in the state vector.
    %
    %   il (r/w, default: [])
    %       Vector of the indices to the linear states in the state vector.
    %
    % METHODS
    %   MixedCLGSSModel(fn, An, fl, Al, Q, h, C, R, m0, P0, in, il)
    %       Constructor that initializes the model when the simple form is 
    %       used. In this case, all arguments must be provided.
    %
    %       fn, An, fl, Al, Q, h, C, R 
    %           Function handles of the form @(xn, t, u) corresponding to
    %           the respective function/matrix as defined in the model
    %           above.
    %
    %       m0, P0
    %           Initial state distribution mean and covariance.
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
    %   AWGNModel, RBGF
    %
    % VERSION
    %   2017-01-02
    % 
    % AUTHORS
    %   Roland Hostettler <roland.hostettler@aalto.fi>
    
    %% Properties
    properties
        % Indices of the linear- and non-linear state components.
        % Preferrably gapless, e.g. in = 1:4, il = 5:12.
        in = [];
        il = [];
    end
    
    %% Private Properties
    % Only to be used in the 'simple' form
    properties (Access = protected)
        % Dynamics
        m_fn = @(xn, t, u) [];
        m_fl = @(xn, t, u) [];
        m_An = @(xn, t, u) [];
        m_Al = @(xn, t, u) [];
        %m_Q = @(xn, t, u) [];
        
        % Observation
        m_h = @(xn, t, u) [];
        m_C = @(xn, t, u) [];
        %m_R = @(xn, t, u) [];
    end

    %% General Methods
    % These should not be overwritten by subclasses
    methods (Access = public)
        %% Constructor
        % Only to be used for 'simple' models
        function self = MixedCLGSSModel(varargin)
            % Input order is: fn, An, fl, Al, Q, h, C, R, m0, P0, in, il
            if nargin == 12
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
            elseif nargin ~= 0
                error('Model not fully specified.');
            end
        end

        %% Process Dynamics
        function x_p = f(self, x, q, t, u)
            xn = x(self.in);
            xl = x(self.il);
            
            xn_p = self.fn(xn, t, u) + self.An(xn, t, u)*xl + q(self.in);
            xl_p = self.fl(xn, t, u) + self.Al(xn, t, u)*xl + q(self.il);
            
            x_p = zeros(size(x));
            x_p(self.in) = xn_p;
            x_p(self.il) = xl_p;
        end
        
        %% Process Noise Covariance
        function Q = Q(self, xn, t, u)
            Q = zeros(length(self.in)+length(self.il));
            Q(self.in, self.in) = self.Qn(xn, t, u);
            Q(self.il, self.il) = self.Ql(xn, t, u);
            Q(self.in, self.il) = self.Qnl(xn, t, u);
            Q(self.il, self.in) = self.Qnl(xn, t, u).';
        end
        
        %% Measurement Function
        function y = g(self, x, r, t, u)
            xn = x(self.in);
            xl = x(self.il);
            y = self.h(xn, t, u) + self.C(xn, t, u)*xl + r;
        end        
    end
    
    %% Default Methods
    % These are only provided for simple models; must be overwritten when
    % subclassing.
    methods (Access = public)
        function fn = fn(self, x, t, u)
            fn = self.m_fn(x, t, u);
        end
        
        function An = An(self, x, t, u)
            An = self.m_An(x, t, u);
        end
        
        function fl = fl(self, x, t, u)
            fl = self.m_fl(x, t, u);
        end
        
        function Al = Al(self, x, t, u)
            Al = self.m_Al(x, t, u);
        end
        
        function Qn = Qn(self, x, t, u)
            Q = self.m_Q(x, t, u);
            Qn = Q(self.in, self.in);
        end
        
        function Ql = Ql(self, x, t, u)
            Q = self.m_Q(x, t, u);
            Ql = Q(self.il, self.il);
        end
        
        function Qnl = Qnl(self, x, t, u)
            Q = self.m_Q(x, t, u);
            Qnl = Q(self.in, self.il);
        end
        
        function h = h(self, x, t, u)
            h = self.m_h(x, t, u);
        end
        
        function C = C(self, x, t, u)
            C = self.m_C(x, t, u);
        end
        
        function R = R(self, x, t, u)
            R = self.m_R(x, t, u);
        end
    end
end
