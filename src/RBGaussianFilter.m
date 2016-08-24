classdef (Abstract) RBGaussianFilter < handle
    % Rao-Blackwellized Gaussian Filter
    % 
    % DESCRIPTION
    %   
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
        model;
        
        % Non-linear state
        m_n;
        P_n;
        
        % Linear state
        m_l;
        P_l;
        
        % Cross-covariance
        Pnl;
        
        % TODO: alias properties for m and n
        
    end

    % Methods
    methods
        function self = RBGaussianFilter(model)
            if nargin == 1
                self.model = model;
            end
        end
        
        function update(self, u, y, t)
            self.timeUpdate(u, t);
            
            % TODO: switch between posterior linearization and regular
            
            self.measurementUpdate(u, y, t);
            self.model.t = t;
        end
        
        function timeUpdate(self, u, t)
            
        end
        
        function measurementUpdate(self, u, y, t)
            
        end
    end
end
