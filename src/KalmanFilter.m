classdef KalmanFilter < handle
    % A generic Unscented Kalman Filter
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
    %   2016-03-02
    % 
    % AUTHORS
    %   Roland Hostettler <roland.hostettler@aalto.fi>   
    
    % Properties
    properties
        % The model of the type 'LGSSModel'
        model;
        
        % Mean and covariance
        m;
        P;
        
        % Predicted mean and covariance
        m_p;
        P_p;
    end

    % Methods
    methods
        %% Constructor
        function self = KalmanFilter(model)
            % Store the model
            self.model = model;
            
            % Initialize the filter
            self.m = model.m0;
            self.P = model.P0;
        end
        
        %% Update the filter
        function update(self, u, y, t)
            self.timeUpdate(u, t);
            self.measurementUpdate(u, y, t);
        end
        
        %% Time update function
        function timeUpdate(self, u, t)
            % Predict the state
            self.m_p = self.model.F*self.m;
            self.P_p = self.model.F*self.P*self.model.F' + self.model.Q;
        end
        
        %% Measurement update function
        function measurementUpdate(self, u, y, t)
                % Kalman gain
                S = self.model.G*self.P_p*self.model.G' + self.model.R;
                K = (self.P_p*self.model.G')/S;
                
                % Time update
                self.m = self.m_p + K*(y - self.model.G*self.m_p);
                self.P = self.P_p - K*S*K';
        end
    end
end
