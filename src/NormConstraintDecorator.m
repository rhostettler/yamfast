classdef NormConstraintDecorator < FilterDecorator
    % Norm-constraints decorator
    % 
    % DESCRIPTION
    %   A decorator implementing norm-constraints on parts of the state.
    %
    % PROPERTIES
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
    %   timeUpdateHook(filter)
    %       Re-implementation of the FilterDecorator's time update hook.
    %
    %   measurementUpdateHook(filter)
    %       Re-implementation of the FilterDecorator's measurement update
    %       hook.
    %
    % SEE ALSO
    %   FilterDecorator
    %
    % VERSION
    %   2016-12-22
    % 
    % AUTHORS
    %   Roland Hostettler <roland.hostettler@aalto.fi>
    
    %% Properties
    properties (Access = public)
        % Constraints
        constraints = {};        
    end

    %% Methods
    methods (Access = public)        
        %% Prediction Hook
        function timeUpdateHook(self, filter)
            filter.m_p = self.normalize(filter.m_p);
        end
        
        %% Measurement Update Hook
        function measurementUpdateHook(self, filter)
            filter.m = self.normalize(filter.m);
        end
    end
    
    %% Private Methods
    methods (Access = private)
        %% Function to Impose the Norm-Constraints
        function m = normalize(self, m)
            Nc = size(self.constraints, 1);
            for n = 1:Nc
                indices = self.constraints{n, 1};
                l = self.constraints{n, 2};
                m(indices) = l/norm(m(indices))*m(indices);
            end
        end
    end
end
