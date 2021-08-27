classdef ExtendedKF < GaussianFilter
    % Extended Kalman filter
    % 
    % DESCRIPTION
    %   Extended Kalman filter where the non-linear state transition
    %   function f and measurement function g are linearized around m_p and
    %   m, respectively, for updating the prediction and posterior
    %   covariances. See [1] for details.
    %
    % PROPERTIES
    %   Does not provide properties of its own, see the GaussianFilter's
    %   documentation for inherited properties.
    %
    % METHODS
    %   Implements all the function as required (and described) by the
    %   GaussianFilter-class. See its documentation for details.
    %
    %   ExtendedKF(model)
    %       Constructor.
    %
    %       model (optional)
    %           State-space model for this filter to operate on. While
    %           specifying a model is optional, omitting the model is only
    %           meaningful in cases where this filter is used together with
    %           an IMM (see IMMFilter). In all other cases, the model must
    %           be specified.
    %
    % REFERENCES
    %   [1] S. S?rkk?, "Bayesian Filtering and Smoothing", Cambridge
    %       University Press, 2013.
    %
    % SEE ALSO
    %   GaussianFilter, KalmanFilter, GenericGF, AdditiveGF
    %
    % VERSION
    %   2016-12-22
    % 
    % AUTHORS
    %   Roland Hostettler <roland.hostettler@aalto.fi>   
    
    %#ok<*PROPLC>
    
    %% Properties
    properties (Access = public)
        % No own properties.
    end

    %% Methods
    methods (Access = public)
        %% Constructor
        function self = ExtendedKF(model)
            if nargin >= 1
                self.model = model;
                self.initialize();
            end
        end
                
        %% Time update function
        function [m_p, P_p] = timeUpdate(self, t, u)
            Q = self.model.Q(self.m, t, u);
            Nq = size(Q, 1);

            % Propagate the state and get the Jacobians w.r.t. x (Fx) and q
            % (Fr)
            [m_p, Fx, Fr] = self.model.f(self.m, zeros(Nq, 1), t, u);
            
            % Calculate the predicted covariance
            P_p = Fx*self.P*Fx' + Fr*Q*Fr';
            P_p = (P_p + P_p')/2;
            
            % Store
            self.m_p = m_p;
            self.P_p = P_p;
        end
        
        %% Measurement update function
        function [m, P] = measurementUpdate(self, y, t, u)
            R = self.model.R(self.m_p, t, u);
            Nr = size(R, 1);

            % Predict the output and get the Jacobians w.r.t. x (Gx) and r
            % (Gr)
            [y_p, Gx, Gr] = self.model.g(self.m_p, zeros(Nr, 1), t, u);

            % Kalman Gain
            S = Gx*self.P_p*Gx' + Gr*R*Gr';
            K = (self.P_p*Gx')/S;

            % Update
            v = y - y_p;
            m = self.m_p + K*v;
            P = self.P_p - K*S*K';
            P = (P + P')/2;
            
            % Store
            self.m = m;
            self.P = P;
            self.y_p = y_p;
            self.S = S;
        end
    end
end
