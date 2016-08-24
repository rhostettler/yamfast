classdef LGSSModel < GenericStateSpaceModel
    % Linear, Gaussian State-Space Model
    % 
    % DESCRIPTION
    %   Implements the state transition and observation functions for the
    %   linear, Gaussian state-space model of the form
    %
    %       x(k) = F(t(k)) x(k-1) + B(t(k)) u(k-1) q(k-1)
    %       y(k) = G(t(k)) x(k) + D(t(k)) u(k) + r(k)
    %
    %   where:
    %       * q(k) ~ N(0, Q)
    %       * r(k) ~ N(0, R)
    %       * x(0) ~ N(m0, P0)
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
    %   2016-03-02
    % 
    % AUTHORS
    %   Roland Hostettler <roland.hostettler@aalto.fi>   
    
    % Properties
    properties
        % System matrices
        F;
        G;
        B = 0;
        D = 0;
    end

    % Methods
    methods
        function self = LGSSModel(F, G, Q, R, m0, P0)
            self.F = F;
            self.G = G;
            self.Q = Q;
            self.R = R;
            self.m0 = m0;
            self.P0 = P0;
        end
        
        function x_p = f(self, x, u, q, t)
            x_p = self.F*x + self.B*u + q;
        end
        
        function y_p = g(self, x, u, r, t)
            y_p = self.G*x + self.D*u + r;
        end
    end
end
