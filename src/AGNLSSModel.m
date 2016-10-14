classdef (Abstract) AGNLSSModel < handle
    % Gaussian Non-Linear State-Space Model
    % 
    % DESCRIPTION
    %   
    %
    %       x[t] = f(x[t-1], u[t], t) + q[t]
    %       y[t] = g(x[t], u[t], t) + r[t]
    %       x[0] ~ N(mu0, P0)
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
    %   2016-
    % 
    % AUTHORS
    %   Roland Hostettler <roland.hostettler@aalto.fi>   
    
    % Properties
    properties
        
        mu0;
        P0;
        
        Q;
        R;
    end

    % Methods
    methods (Abstract)
        f(self, x, u, q, t);
        g(self, x, u, r, t);        
    end
        
    methods (Access = public)
        function x0 = px0_rand(self, M)
            x0 = self.mu0*ones(1, M) ...
                + chol(self.P0).'*randn(size(self.P0, 1), M);
        end
        
        function xp = px_rand(self, x, u, t)
            q = chol(self.Q).'*randn(size(self.Q, 1), size(x, 2));
            xp = self.f(x, u, q, t);
        end
        
        function py = py_eval(self, y, x, u, t)
            py = mvnpdf(y.', self.g(x, u, t).', self.R).';
        end
    end
end
