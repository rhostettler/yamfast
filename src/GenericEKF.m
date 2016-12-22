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
            
            warning('The interface of this filter has changed, make sure you update your code.');
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
        function update(self, y, t, u)
            self.timeUpdate(t, u);
            self.measurementUpdate(y, t, u);
        end
        
        %% Time update function
        function timeUpdate(self, t, u)
            Q = self.model.Q(self.m, t, u);
            Nv = size(Q, 1);

            % Propagate the state and get the Jacobians w.r.t. x(t-1) (F) 
            % and v(t) (B)
            [m_p, Fx, Fr] = self.model.f(self.m, zeros(Nv, 1), t, u);
            
            % Calculate the predicted covariance
            P_p = Fx*self.P*Fx' + Fr*Q*Fr';
                        
            % Store
            self.m_p = m_p;
            self.P_p = (P_p + P_p')/2;
        end
        
        %% Measurement update function
        function measurementUpdate(self, y, t, u)
            R = self.model.R(self.m_p, t, u);
            Nw = size(R, 1);

            % Predict the output and get the Jacobians w.r.t. x(t) (G) and
            % w(t) (D)
            [y_p, Gx, Gr] = self.model.g(self.m_p, zeros(Nw, 1), t, u);

            % Kalman Gain
            S = Gx*self.P_p*Gx' + Gr*R*Gr';
            K = (self.P_p*Gx')/S;

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
