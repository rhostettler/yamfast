classdef (Abstract) CLGSSModel < GenericStateSpaceModel
    % Mixed linear/non-linear Gaussian State Space Model Template Class
    % 
    % DESCRIPTION
    %   A generic mixed linear/non-linear Gaussian state space model of the
    %   form
    %
    %       xn[t] = fn(xn[t-1]) + An(xn[t-1]) xl[t-1] + qn[t]
    %       xl[t] = fl(xn[t-1]) + Al(xn[t-1]) xl[t-1] + ql[t]
    %       y[t] = h(xn[t]) + C(xn[t]) xl[t-1] + r[t]
    %
    %   where q[t] = [qn[t], ql[t]]^T ~ N(0, Q(xn[t-1]), r[t] ~ N(0,
    %   R(xn[t])), and x[0] = [xn[0], xl[0]]^T ~ N(m0, P0).
    %
    % PROPERTIES
    %   
    %
    % METHODS
    %
    %
    % SEE ALSO
    %   GenericCLGSSModel, RBGF
    %
    % VERSION
    %   2016-10-14
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

    %% Abstract Methods
    % These need to be implemented by subclasses
    methods (Abstract)        
        Qn(self, xn, u, t);
        Ql(self, xn, u, t);
        Qnl(self, xn, u, t);
        fn(self, xn, u, t);
        An(self, xn, u, t);
        fl(self, xn, u, t);
        Al(self, xn, u, t);
        h(self, xn, u, t);
        C(self, xn, u, t);
        R(self, xn, u, t);
    end
    
    %% General Methods
    methods (Access = public)  
        function x_p = f(self, x, u, q, t)
            xn = x(self.in);
            xl = x(self.il);
            
            xn_p = self.fn(xn, u, t) + self.An(xn, u, t)*xl + q(self.in);
            xl_p = self.fl(xn, u, t) + self.Al(xn, u, t)*xl + q(self.il);
            
            x_p = zeros(size(x));
            x_p(self.in) = xn_p;
            x_p(self.il) = xl_p;
        end
        
        function y = g(self, x, u, r, t)
            xn = x(self.in);
            xl = x(self.il);
            y = self.h(xn, u, t) + self.C(xn, u, t)*xl + r;
        end
    end
end
