classdef RTSSmoother < KalmanFilter
    % Rauch-Tung-Striebel smoother
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
    %   2017-01-05
    % 
    % AUTHORS
    %   Roland Hostettler <roland.hostettler@aalto.fi>   
    
    % TODO:
    %   * I should define a generic smoother interface that does all the
    %     stuff in the smooth()-function ('GaussianSmoother, subclass of
    %     'GaussianFilter'?), including possible decorators
    
    % Properties
    properties
        % Filtered mean and covariance
        m_f;
        P_f;
        
        % Predicted mean and covariance
        mp_f;
        Pp_f;
        
        % Smoothed mean and covaraiance
        m_s;
        P_s;
    end

    % Methods
    methods
        %% Constructor
        function self = RTSSmoother(model)
            self = self@KalmanFilter(model);
        end
        
        %% 
        function [m_s, P_s] = smooth(self, y, t, u)
            Nx = size(self.m, 1);
            N = size(y, 2);
            
            self.m_f = zeros(Nx, N);
            self.P_f = zeros(Nx, Nx, N);
            self.mp_f = zeros(Nx, N);
            self.Pp_f = zeros(Nx, Nx, N);
            self.m_s = zeros(Nx, N);
            self.P_s = zeros(Nx, Nx, N);
            
            %% Sanity Checks
            if isempty(u)
                u = zeros(1, N);
            end

            %% Forward-Filter
            for n = 1:N
                self.forwardIteration(y(:, n), t(n), u(:, n), n);
            end
            
            %% Backward-Filter & Smoothing
            self.initializeBackwardIteration(y(:, N), t(N), u(:, N));
            for n = N-1:-1:1
                self.backwardIteration(y(:, n), t(n), u(:, n), n);
            end
            
            m_s = self.m_s;
            P_s = self.P_s;
        end    
    end
    
    methods (Access = private)
        function forwardIteration(self, y, t, u, n)
            self.update(y, t, u);
            self.m_f(:, n) = self.m;
            self.P_f(:, :, n) = self.P;
            self.mp_f(:, n) = self.m_p;
            self.Pp_f(:, :, n) = self.P_p;
        end
        
        function initializeBackwardIteration(self, y, t, u)
            self.m_s(:, end) = self.m_f(:, end);
            self.P_s(:, :, end) = self.P_f(:, :, end);
        end
        
        function backwardIteration(self, y, t, u, n)
            P = self.P_f(:, :, n);
            Pp = self.Pp_f(:, :, n+1);
            F = self.model.F(t, u);
            
            L = (P*F')/Pp;
            self.m_s(:, n) = ( ...
                self.m_f(:, n) + L*(self.m_s(:, n+1) - self.mp_f(:, n+1)) ...
            );
            self.P_s(:, :, n) = P + L*(self.P_s(:, :, n+1) - Pp)*L';
        end
    end
end
