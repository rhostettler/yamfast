classdef (Abstract) GaussianFilter < handle
    % General Gaussian filter
    % 
    % DESCRIPTION
    %   Interface definition for Gaussian (assumed density) filters such as
    %   the Kalman filter (KF), extended KF (EKF), or unscented Kalman
    %   filter (UKF).
    %
    %   Gaussian filters approximate the prediction- and filtering
    %   densities as
    %
    %       p(x[n] | y_[1:n-1]) = N(x[n]; m_p, P_p)
    %       p(x[n] | y_[1:n]) = N(x[n]; m, P)
    %
    %   (which is of course exact for the Kalman filter).
    %
    %   This class provides all the properties and functionality common to
    %   all Gaussian filters.
    %
    % PROPERTIES
    %   model (r/w)
    %       The model used in this filter.
    %
    %   m, P, D (r/w)
    %       Posterior mean, covariance, and cross-covariacne from the last 
    %       measurement update.
    %
    %   m_p, P_p, C (r/w)
    %       Predicted mean, covariance, and cross-covariance from the last 
    %       time update.
    %
    %   y_p, S (r/w)
    %       Innovation and innovation covariance.
    %
    % METHODS
    %   [m, P] = update(y, t, u)
    %       Complete filter update including time- and measurement update.
    %       Returns the estimated mean and covariance (which can also be
    %       accessed through the properties 'm' and 'P' after the update.
    %
    %       y   Measurement vector (Ny x 1).
    %       t   Time of measurement.
    %       u   Control input.
    %
    %   [m_p, P_p] = timeUpdate(t, u)
    %       Time update step (prediction). This is a virtual function and
    %       must be implemented by the actual filtering algorithm.
    %
    %       t   Time.
    %       u   Control input.
    %
    %   [m, P] = measurementUpdate(y, t, u)
    %       Measurement update step. This is a virtual function and must be
    %       implemented by the actual filtering algorithm.
    %
    %       y   Measurement vector (Ny x 1).
    %       t   Time of measurement.
    %       u   Control input.
    %
    %   initialize()
    %       (Re-)initializes the filter.
    %
    % REFERENCES
    %   [1] S. S??rkk??, "Bayesian Filtering and Smoothing", Cambridge
    %       University Press, 2013.
    %
    % SEE ALSO
    %   KalmanFilter, ExtendedKF, GenericGF, AdditiveGF
    %
    % VERSION
    %   2016-12-22
    % 
    % AUTHORS
    %   Roland Hostettler <roland.hostettler@aalto.fi>   
    
    %#ok<*PROPLC>
    
    %% Properties
    properties (Access = public)
        % The model
        model = [];
        
        % Mean and covariance
        m;
        P;
        D;
        
        % Predicted mean and covariance
        m_p;
        P_p;
        C;
        
        % Innovation and innovation covariance
        y_p;
        S;
    end
    
    %% Protected Properties
    properties (Access = protected)
        decorators;
    end

    %% Methods (Prototypes)
    % These are implemented by the actual filters
    methods
        timeUpdate(self, t, u);
        measurementUpdate(self, y, t, u);
    end
    
    %% Methods
    methods (Access = public)
        %% Update the filter
        function [m, P] = update(self, y, t, u)
            self.timeUpdate(t, u);
            self.timeUpdateHook();
            self.measurementUpdate(y, t, u);
            self.measurementUpdateHook();
            m = self.m;
            P = self.P;
        end
        
        %% Run Update Hooks
        function timeUpdateHook(self)
            for i = 1:length(self.decorators)
                self.decorators{i}.timeUpdateHook(self);
            end
        end
        
        function measurementUpdateHook(self)
            for i = 1:length(self.decorators)
                self.decorators{i}.measurementUpdateHook(self);
            end
        end
        
        %% Filter a Batch of Data
        function [m, P] = filter(self, y, t, u)
            N = size(y, 2);
            if nargin < 4
                u = zeros(1, N);
            end
            
            model = self.model;
            Nx = size(model.m0, 1);
            m = zeros(Nx, N);
            P = zeros(Nx, Nx, N);
            
            for n = 1:N
                self.update(y(:, n), t(n), u(n));
                m(:, n) = self.m;
                P(:, :, n) = self.P;
                model.t = t(n);
            end            
%             self.m = m;
%             self.P = P;
        end
        
        %% Adds a Decorators
        function addDecorator(self, decorator)
            self.decorators{length(self.decorators)+1} = decorator;
        end
        
        %% (Re-)Initialize the Filter
        function initialize(self)
            if ~isempty(self.model)
                self.m = self.model.m0;
                self.P = self.model.P0;
            end
        end
    end
end
