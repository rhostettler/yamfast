classdef NormConstrainedURTSS < NormConstrainedUKF
    % Norm-constrained Unscented Rauch-Tung-Striebel Smoother
    % 
    % DESCRIPTION
    %   Unscented Rauch-Tung-Striebel smoother with brute-force
    %   normalization. Note that the brute-force normalization is probabliy
    %   not optimal in any sense.
    %
    % PROPERTIES
    %   alpha (default: 1e-3)
    %       'alpha' filter parameter.
    %   
    %   beta (default: 2)
    %       'beta' filter parameter.
    %
    %   kappa (default: 0)
    %       'kappa' filter parameter.
    %
    %   mf, Pf
    %       Filter mean and covariance.
    %
    %   mpf, Ppf
    %       Predicted mean and covaraince.
    %
    %   ms, Ps
    %       Smoothed mean and covariance.
    %
    %   constraints
    %       An array with the norm constraints. Each row in the array is
    %       one constraint consisting of the affected state indices and the
    %       constraint value, for example, if
    %
    %           constraints = {1:3, 4}
    %
    %       then
    %
    %           x(1:3) = 4/norm(x(1:3))*x(1:3);
    %
    %   model
    %       The model associated with this smoother.
    %
    % METHODS
    %   smooth(u, y, t)
    %       Smooth a batch of data with
    %
    %           u   Control inputs, Nu x N
    %           y   Measurements, Ny x N
    %           t   Time, 1 x N
    %
    % SEE ALSO
    %   NormConstrainedUKF
    %
    % VERSION
    %   2016-06-01
    % 
    % AUTHORS
    %   Roland Hostettler <roland.hostettler@aalto.fi>
    
    % Disable some mlint messages
    %#ok<*PROPLC>
    
    %% Properties
    properties
        % Mean & covariance
        mf;
        Pf;
        
        mpf;
        Ppf;
        
        % Smoothed mean & covariance
        ms;
        Ps;                
    end

    %% Methods
    methods (Access = public)
        %% Constructor
        function self = NormConstrainedURTSS(varargin)
            if nargin >= 1
                self.model = varargin{1};
                self.m = self.model.m0;
                self.P = self.model.P0;
            end
        end
        
        %% Filter Update
        % Input
        %   u: Control signal
        %   y: Measurement
        %   t: Time
        function smooth(self, u, y, t)
            N = length(t);
            Nx = size(self.model.m0, 1);
            Nq = size(self.model.Q([], [], 0), 1);
            
            % Preallocate
            self.mf = zeros(Nx, N);
            self.Pf = zeros(Nx, Nx, N);
            self.mpf = zeros(Nx, N);
            self.Ppf = zeros(Nx, Nx, N);
            self.Cp = zeros(Nx+Nq, Nx+Nq, N);
            
            %% Filtering Pass
            for n = 1:N
                self.update(u(:, n), y(:, n), t(n));
                self.model.t = t(n);
                
                % store to local
                self.mf(:, n) = self.m;
                self.Pf(:, :, n) = self.P;
                self.mpf(:, n) = self.mp;
                self.Ppf(:, :, n) = self.Pp;
                self.Cp(:, :, n) = self.Cp;
            end
            
            %% Initialize the Smoother
            ms = self.m;
            Ps = self.P;
            
            self.ms(:, N) = ms;
            self.Ps(:, :, N) = Ps;
            
            %% Smoothing Pass
            for n = N-1:-1:1
                % Smoothing
                C = self.Cp(1:Nx, :, n+1);
                D = C/self.Ppf(:, :, n+1);
                ms = self.mf(:, n) + D*(ms - self.mpf(:, n+1));
                Ps = self.Pf(:, :, n) + D*(Ps - self.Ppf(:, :, n+1))*D';                
                
                % Normalize
                ms = self.normalize(ms);
                
                % Store
                self.Ps(:, :, n) = Ps;
                self.ms(:, n) = ms;
            end
        end
    end
end
