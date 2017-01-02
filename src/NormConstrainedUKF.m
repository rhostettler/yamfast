classdef NormConstrainedUKF < GenericUKF
    % Norm-constrained unscented Kalman filter
    % 
    % DESCRIPTION
    %   Unscented Kalman filter where part of the state is norm
    %   constrained.
    %
    % PROPERTIES
    %   Inherits all the properties of GaussianFilter and has the following
    %   properties of its own.
    %
    %   alpha (r/w, default: 1)
    %       The UT's 'alpha' parameter.
    %   
    %   beta (r/w, default: 0)
    %       The UT's 'beta' parameter.
    %
    %   kappa (r/w, default: 1)
    %       The UT's 'kappa' parameter.
    %
    %   constraints
    %       An array with the norm constraints. Each row in the array is
    %       one constraint consisting of the affected state indices and the
    %       constraint value, for example, if
    %
    %           constraints = {1:3, 4}
    %
    %       then
    %
    %           x(1:3) = 4/norm(x(1:3))*x(1:3);
    %
    % METHODS
    %   Inherits all of GaussianFilter's methods.
    %
    % SEE ALSO
    %   GaussianFilter, GenericGF, UnscentedTransform, 
    %   NormConstraintDecorator, NormConstrainedURTSS
    %
    % VERSION
    %   2017-01-02
    % 
    % AUTHORS
    %   Roland Hostettler <roland.hostettler@aalto.fi>
        
    %% Mapped Properties
    properties (Dependent)
        % Constraints
        constraints;
    end
    
    %% Private Properties
    properties (Access = private)
        constraintDecorator = NormConstraintDecorator();
    end
    
    %% Property Access Methods
    methods
        function set.constraints(self, constraints)
            self.constraintDecorator.constraints = constraints;
        end
        
        function constraints = get.constraints(self)
            constraints = self.constriantDecorator.constraints;
        end
    end

    %% Methods
    methods (Access = public)
        %% Constructor
        function self = NormConstrainedUKF(model)
            if nargin < 1
                model = [];
            end            
            self@GenericUKF(model);
            self.addDecorator(self.constraintDecorator);
        end
    end
end
