classdef GenericUKF < GenericGF
    % A generic unscented Kalman filter
    % 
    % DESCRIPTION
    %   A generic unscented Kalman filter for models of the form
    %
    %       x[n] = f(x[n-1], q[n], t[n], u[n])
    %       y[n] = g(x[n], r[n], t[n], u[n])
    %
    %   with x[0] ~ N(m[0], P[0]), q[n] ~ N(0, Q[n]), r[n] ~ N(0, R[n]).
    %
    % PROPERTIES
    %   Inherits the properties provided by GaussianFilter as well as the
    %   following additional properties.
    %
    %   alpha (r/w, default: 1)
    %       The unscented transform's 'alpha' parameter.
    %
    %   beta (r/w, default: 1)
    %       The unscented transform's 'beta' parameter.
    %
    %   kappa (r/q, default: 1)
    %       The unscented transform's 'kappa' parameter.
    %
    % METHODS
    %   Inherits the methods provided by GaussianFilter, see its
    %   documentation for details.
    %
    % SEE ALSO
    %   GaussianFilter, GenericGF, GenericEKF
    %
    % VERSION
    %   2017-01-02
    % 
    % AUTHORS
    %   Roland Hostettler <roland.hostettler@aalto.fi>
    
    %% Properties
    properties (Dependent)
        alpha;
        beta;
        kappa;
    end
    
    %% Property Access Methods
    methods
        function set.alpha(self, alpha)
            self.rule.alpha = alpha;
        end
        
        function alpha = get.alpha(self)
            alpha = self.rule.alpha;
        end
        
        function set.beta(self, beta)
            self.rule.beta = beta;
        end
        
        function beta = get.beta(self)
            beta = self.rule.beta;
        end
        
        function set.kappa(self, kappa)
            self.rule.kappa = kappa;
        end
        
        function kappa = get.kappa(self)
            kappa = self.rule.kappa;
        end    
    end

    %% Methods
    methods (Access = public)
        %% Constructor
        function self = GenericUKF(model)
            if nargin < 1
                model = [];
            end
            self@GenericGF(model);
        end
    end
end
