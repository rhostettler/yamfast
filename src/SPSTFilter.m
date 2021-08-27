classdef SPSTFilter < handle
    % Sigma-Point Student's t Filter
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
        model;
        
        % Mean and covariance
        m;
        P;
        v = 1e3;  % Set to "inf" because, why not?
        
        % Predicted mean and covariance
        m_p;
        P_p;
        v_p;
        
        % Weight parameters
        kappa = 0;
        d = 3;
        
        % Needed to calculate the likelihood
        Pyy
    end

    %% Methods
    methods (Access = public)
        %% Constructor
        function self = SPSTFilter(model)
            % 
            if nargin >= 1
                % Store the model
                self.model = model;

                % Initialize the filter
                self.m = model.m0;
                self.P = model.P0;
                self.v = model.v0;
            end
        end
        
        %% Update the filter
        function update(self, u, y, t)
            self.timeUpdate(u, t);
            self.measurementUpdate(u, y, t);
        end
        
        %% Time update function
        function timeUpdate(self, u, t)
           % TODO: Make this generic/selectable
            v = min([self.v, self.model.vq, self.model.vr]);
            
            Q = (v - 2)/v*self.model.Q;
            P = (v - 2)*self.v/v/(self.v - 2)*self.P;
            P = (P+P')/2;
            f = @(x, n) self.model.f(x, u, n, t);
            
            % Predict
            [Wxn, Sxn] = self.weights(self.d, size(self.m, 1)+size(Q, 1), v);
            [m_p, P_p, v_p] = stf_predict(self.m, P, v, f, Q, Sxn, Wxn);
            
            % Store
            self.m_p = m_p;
            self.P_p = (P_p+P_p')/2;
            self.v_p = v_p;
        end
        
        %% Measurement update function
        function measurementUpdate(self, u, y, t)            
           % TODO: Make this generic/selectable
            v = min([self.v, self.model.vq, self.model.vr]);
            R = (v-2)/v*self.model.R;
            h = @(x, e) self.model.g(x, u, e, t);
                
            % Update
            [Wxe, Sxe] = self.weights(self.d, size(self.m, 1)+size(R, 1), v);
            [m, P, v, Pyy] = stf_update(y, self.m_p, self.P_p, self.v_p, h, R, Sxe, Wxe);
            
            % Store
            self.m = m;
            self.P = (P+P')/2;
            self.v = v;
            self.Pyy = Pyy;
        end
       
        %% Calculates the sigma points and their weights
        % d: Degrees of integration
        % L: Augmented state dimension
        % v: Degrees of freedom
        function [W, X] = weights(self, d, L, v)
            kappa = self.kappa;
            switch d
                case 3
                    [W, X] = st3_ws(L, v, kappa);
                case 5
                    [W, X] = st5_ws(L, v, kappa);
                case 7
                    [W, X] = st7_ws(L, v, kappa);
                otherwise
                    error('Degree of integration must be 3, 5 or 7.')
            end
        end
    end
end
