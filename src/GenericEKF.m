classdef GenericEKF < handle
    % Generic extended Kalman filter
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
        % The model of the type 'GenericStateSpaceModel'
        model = [];
        
        % Mean and covariance
        m;
        P;
        
        % Predicted mean and covariance
        m_p;
        P_p;
        
        % Innovation and innovation covariance
        v;
        Pyy;
    end

    % Methods
    methods
        %% Constructor
        function self = GenericEKF(model)
            
            if nargin >= 1
                % Store the model
                self.model = model;

                % Initialize the filter
                self.reset();
%                 self.m = model.m0;
%                 self.P = model.P0;
            end
        end
        
        %% (Re-)Initialize
        function reset(self)
            if ~isempty(self.model)
                self.model.t = 0;
                self.m = self.model.m0;
                self.P = self.model.P0;
            end
        end
        
        %% Update the filter
        function update(self, u, y, t)
            self.timeUpdate(u, t);
            self.measurementUpdate(u, y, t);
            self.model.t = t;
        end
        
        %% Time update function
        function timeUpdate(self, u, t)
            Q = self.model.Q;
            Nv = size(Q, 1);

            % Propagate the state and get the Jacobians w.r.t. x(t-1) (F) 
            % and v(t) (B)
            [m_p, F, B] = self.model.f(self.m, u, zeros(Nv, 1), t);
            
            % Calculate the predicted covariance
            P_p = F*self.P*F' + B*Q*B';
                        
            % Store
            self.m_p = m_p;
            self.P_p = (P_p + P_p')/2;
        end
        
        %% Measurement update function
        function measurementUpdate(self, u, y, t)
            R = self.model.R;
            Nw = size(R, 1);

            % Predict the output and get the Jacobians w.r.t. x(t) (G) and
            % w(t) (D)
            [y_p, G, D] = self.model.g(self.m_p, u, zeros(Nw, 1), t);

            % Kalman Gain
            S = G*self.P_p*G' + D*R*D';
            K = (self.P_p*G')/S;

            % Update
            v = y - y_p;
            m = self.m_p + K*v;
            P = self.P_p - K*S*K';

            % Store
            self.m = m;
            self.P = (P + P')/2;
            self.v = v;
            self.Pyy = S;
        end
    end
end
