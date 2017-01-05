classdef KalmanFilter < GaussianFilter
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
    end

    % Methods
    methods
        %% Constructor
        function self = KalmanFilter(model)
            if nargin == 1
                self.model = model;
                self.initialize();
            end
        end
                
        %% Time update function
        function [m_p, P_p] = timeUpdate(self, t, u)
            m_p = self.model.F(t, u)*self.m;
            P_p = ( ...
                self.model.F(t, u)*self.P*self.model.F(t, u)' ...
                + self.model.Q([], t, u) ...
            );
            
            self.m_p = m_p;
            self.P_p = P_p;
        end
        
        %% Measurement update function
        function [m, P] = measurementUpdate(self, y, t, u)
                m_p = self.m_p;
                P_p = self.P_p;
                G = self.model.G(t, u);
                
                S = G*P_p*G' + self.model.R([], t, u);
                K = (P_p*G')/S;
                v = (y - G*m_p);
                m = m_p + K*v;
                P = P_p - K*S*K';
                
                self.m = m;
                self.P = P;
                self.v = v;
                self.Pyy = S;
        end
    end
end
