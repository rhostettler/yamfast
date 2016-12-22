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
    %   m, P (r/w)
    %       Posterior mean and covariance from the last measurement update.
    %
    %   m_p, P_p (r/w)
    %       Predicted mean and covariance from the last time update.
    %
    %   v, Pyy (r/w)
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
    % REFERENCES
    %   [1] S. S?rkk?, "Bayesian Filtering and Smoothing", Cambridge
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
        
        % Predicted mean and covariance
        m_p;
        P_p;
        
        % Innovation and innovation covariance
        v;
        Pyy;
    end

    %% Methods (Prototypes)
    % These are implemented by the actual filters
    methods
        timeUpdate(t, u);
        measurementUpdate(y, t, u);
    end
    
    %% Methods
    methods (Access = public)
        %% Update the filter
        function [m, P] = update(self, y, t, u)
            self.timeUpdate(t, u);
            [m, P] = self.measurementUpdate(y, t, u);
        end
    end
    
    %% Helper Methods
    methods (Access = protected)
        %% (Re-)Initialize
        function initialize(self)
            if ~isempty(self.model)
                self.m = self.model.m0;
                self.P = self.model.P0;
            end
        end
    end
end
