classdef GaussianSmoother < handle
    % Generic Gaussian smoother
    % 
    % DESCRIPTION
    %   Gaussian assumed density smoother for non-linear state-space
    %   models. Approximates the smoothing density
    %
    %       p(x_n | y_{1:N}) = N(m_s, P_s).
    %
    %   The smoother implements the RTS recursion and can use any Gaussian
    %   filter (e.g. UKF, EKF, etc.) to calculate the smoothed mean and
    %   covariance. By default, an unscented Kalman filter is used.
    %
    % PROPERTIES
    %   m_s, P_s (r/w)
    %       Smoothed mean and covariance.
    %
    %   E (r/w)
    %       Cross-covariance matrix.
    %
    % METHODS
    %   GaussianSmoother(model, filter)
    %       Initializes the smoother.
    %
    %       model
    %           The state-space model for this smoother to operate on.
    %
    %       filter (optional, default: GenericGF)
    %           The filter to use in the forward pass. If omitted, a filter
    %           of the type GenericGF is used which defaults to an UKF,
    %           yielding an URTSS.
    %   
    %   [m_s, P_s] = smooth(y, t, u)
    %       Smooths the batch of data y and returns the smoothed mean and
    %       covariance.
    %
    %       y   Data vector.
    %       t   Time vector.
    %       u   Control input (optional)
    %
    % SEE ALSO
    %   GaussianFilter, GenericGF, AdditiveGF
    %
    % VERSION
    %   2017-01-17
    % 
    % AUTHORS
    %   Roland Hostettler <roland.hostettler@aalto.fi>   
    
    %#ok<*PROPLC>
    
    %% Properties
    properties (Access = public)
        % Smoothed mean, covariance, and cross-covariance
        m_s;
        P_s;
        E;
    end
    
    properties (Access = private)
        filter = [];
        m;
        P;
        m_p;
        P_p;
        C;
    end

    %% Public Methods
    methods (Access = public)
        %% Initialization
        function self = GaussianSmoother(model, filter)
            if nargin < 2
                filter = GenericGF(model);
            else
                filter.model = model;
                filter.initialize();
            end
            self.filter = filter;
        end
        
        %% Smoothing Function
        function [m_s, P_s] = smooth(self, y, t, u)
            filter = self.filter;
            Nx = size(filter.model.m0, 1);
            N = length(t);
            
            % Preallocate
            m = zeros(Nx, N);
            P = zeros(Nx, Nx, N);
            m_p = m;
            P_p = P;
            C = P;
            m_s = m;
            P_s = P;
            E = P;
            
            % Initialize an empty input if not set
            if nargin < 4 || isempty(u)
                u = zeros(1, N);
            end

            %% Forward Pass
            for n = 1:N
                filter.update(y(:, n), t(n), u(:, n));
                m(:, n) = filter.m;
                P(:, :, n) = filter.P;
                m_p(:, n) = filter.m_p;
                P_p(:, :, n) = filter.P_p;
                C(:, :, n) = filter.C;
            end
            
            %% Backward Pass
            % Initialize backward pass
            m_s(:, N) = m(:, N);
            P_s(:, :, N) = P(:, :, N);
            for n = N-1:-1:1
                % TODO: Update Hooks
                G = C(:, :, n+1)/P_p(:, :, n+1);
                m_s(:, n) = m(:, n) + G*(m_s(:, n+1) - m_p(:, n+1));
                P_s(:, :, n) = P(:, :, n) ...
                    + G*(P_s(:, :, n+1) - P_p(:, :, n+1))*G';
                E(:, :, n) = G*P_s(:, :, n+1);
            end
            
            %% Store
            self.m = m;
            self.P = P;
            self.m_p = m_p;
            self.P_p = P_p;
            self.C = C;
            self.m_s = m_s;
            self.P_s = P_s;
            self.E = E;
        end
    end
end
