classdef RBFFBSi < handle
    % Rao-Blackwellized forward-filtering backward-simulation smoother
    % 
    % DESCRIPTION
    %   Rao-Blackwellized forward-filtering backward-simulation particle 
    %   smoother for conditionally linear Gaussian state-space models 
    %   according to [1]
    %
    %   TODO: Enhance documentation
    %
    % PROPERTIES
    % 
    %
    % METHODS
    %
    % REFERENCES
    %   [1] F. Lindsten, P. Bunch, S. Särkkä, T. Schön, S. Godsill,
    %       "Rao-Blackwellized Particle Smoothers for Conditionally Linear
    %       Gaussian Models", Journal of Selected Topics in Signal
    %       Processing, vol. 10, no. 2, March 2016
    %
    % SEE ALSO
    %
    %
    % VERSION
    %   2017-01-17
    % 
    % AUTHORS
    %   Roland Hostettler <roland.hostettler@aalto.fi>   
    
     %#ok<*PROPLC>
     %#ok<*PROP>
    
    %% Properties
    properties (Access = public)
        ss;
        s;
        w;
        z;
        P;
        zs;
        Ps;
    end
    
    %% Aliases
    properties (Dependent)
    end
    
    %% Private Properties
    properties (Access = private)
        M;
        Ms;
        filter;
        model;
        type = 'mixed';
        Omega_hat;
        lambda_hat;
        
        s0;
        
        Omega;
        lambda;
        
        shat_s;
        zhat_s;
        P_s;
    end

    %% Public Methods
    methods (Access = public)
        %% Constructor
        function self = RBFFBSi(model, M, Ms, type)
            % TODO: Make generic
            % TODO: We might want to remove the type later on and determine
            % the type based on the model. (And implement of course).
            switch type
                case 'hierarchical'
                    error('Filter not implemented yet, sorry.');
                case 'mixed'
                    self.filter = MixedRBPF(model, M);
                otherwise
                    error('Unknown model type');
            end
            self.type = type;
            self.M = M;
            self.Ms = Ms;
            self.model = model;
        end
        
        %% Smoothing Function
        function xhat_s = smooth(self, y, t, u)
            %% Preallocate
            M = self.M;
            Ms = self.Ms;
            Ns = size(self.filter.s, 1);
            Nz = size(self.filter.z, 1);
            N = size(y, 2);
            
            % These store the filtered particles
            self.s = zeros(Ns, M, N);
            self.w = zeros(1, M, N);
            self.z = zeros(Nz, M, N);
            self.P = zeros(Nz, Nz, M, N);
            
            % These store the smoothed particles
            self.ss = zeros(Ns, Ms, N);
            self.zs = zeros(Nz, Ms, N);
            self.Ps = zeros(Nz, Nz, Ms, N);
            
            % These store the Omega- and lambda-statistics for the smoothed
            % linear states
            self.lambda = zeros(Nz, Ms, N);
            self.Omega = zeros(Nz, Nz, Ms, N);
            
            % Store the point estimates
            self.shat_s = zeros(Ns, N);
            self.zhat_s = zeros(Nz, N);
            self.P_s = zeros(Nz, Nz, N);
            
            %% Sanity Checks
            if nargin < 4 || isempty(u)
                u = zeros(1, N);
            end

            %% Actual Smoothing
            self.s0 = self.filter.s;
            for n = 1:N
                self.forwardIteration(y(:, n), t(n), u(:, n), n);
            end
            
            % TODO: There's actually something wrong with the timing. We're
            %       extending from n+1 to n, thus, the time should possibly
            %       be t(n+1). Also, we should set the model time to t(n).
            %       It should be similar for u(n).
            self.initializeBackwardIteration(y(:, N), t(N), u(:, N));
            for n = N-1:-1:1
                self.backwardIteration(y(:, n), t(n), u(:, N), n);
            end
            self.resampleInitialParticles();
            self.smoothLinearStates(y, t, u);
            
            %% Point Estimate
            xhat_s = zeros(Ns + Nz, N);
            xhat_s(self.model.in, :) = self.shat_s;
            xhat_s(self.model.il, :) = self.zhat_s;
        end
    end
    
    %% Private Methods
    methods (Access = private)
        %% Forward Iteration
        function forwardIteration(self, y, t, u, n)
            filter = self.filter;
            filter.update(y, t, u);
            self.s(:, :, n) = filter.s;
            self.w(:, :, n) = filter.w;
            self.z(:, :, n) = filter.z;
            self.P(:, :, :, n) = filter.P;
        end
        
        %% Initialization of the Backward Pass
        function initializeBackwardIteration(self, y, t, u)
            Ms = self.Ms;
            Nz = size(self.z, 1);
            lambda_hat = zeros(Nz, Ms);
            Omega_hat = zeros(Nz, Nz, Ms);

            ri = resample(self.w(:, :, end));
            s = self.s(:, ri, end);
            for m = 1:self.Ms
                h = self.model.h(s(:, m), t, u);
                C = self.model.C(s(:, m), t, u);
                R = self.model.R(s(:, m), t, u);
                lambda_hat(:, m) = C'/R*(y - h);
                Omega_hat(:, :, m) = C'/R*C;
            end
            
            self.ss(:, :, end) = s;
            self.shat_s(:, end) = mean(s, 2);
            self.lambda_hat = lambda_hat;
            self.Omega_hat = Omega_hat;
        end
        
        %% Backward Iteration
        function backwardIteration(self, y, t, u, n)
            model = self.model;
            Ms = self.Ms;
            ss = self.ss(:, :, n+1);
            s = self.s(:, :, n);
            w = self.w(:, :, n);
            z = self.z(:, :, n);
            P = self.P(:, :, :, n);
            
            for m = 1:Ms
                % 3a)-3c): Backward prediction
                [wtilde, Omega, lambda] = self.calculateBackwardWeights(ss(:, m), s, w, z, P, t, u);
                
                % 3d)-3f): Extend the state trajectory
                ri = resample(wtilde);
                j = ri(randi(Ms));
                s_n = s(:, j);
                Omega = Omega(:, :, j);
                lambda = lambda(:, j);
                
                % 3g) Compute Omega_hat, lambda_hat
                h = model.h(s_n, t, u);
                C = model.C(s_n, t, u);
                R = model.R(s_n, t, u);
                Omega_hat = Omega + C'/R*C;
                lambda_hat = lambda + C'/R*(y-h);
                
                % Store the extended trajectory and sufficient statistics
                % for this particle
                self.ss(:, m, n) = s_n;
                self.Omega_hat(:, :, m) = Omega_hat;
                self.lambda_hat(:, m) = lambda_hat;
                
                % Store the sufficient statistics
                self.Omega(:, :, m, n) = Omega;
                self.lambda(:, m, n) = lambda;
            end
            
            % Point estimate
            self.shat_s(:, n) = mean(self.ss(:, :, n), 2);
        end

        %% Resamples the Initial Particles
        % Used to initializing smoothing of the linear states in the mixed
        % model
        function resampleInitialParticles(self)
            model = self.model;
            M = self.M;
            Ms = self.Ms;
            s0 = self.s0;
            w0 = 1/M*ones(1, M);
            z0 = model.m0(model.il)*ones(1, M);
            P0 = repmat(model.P0(model.il, model.il), [1, 1, M]);
            
            for m = 1:Ms
                % 3a)-3c): Backward prediction
                wtilde = self.calculateBackwardWeights(self.ss(:, m, 1), s0, w0, z0, P0, [], []);
                
                % 3d)-3f): Extend the state trajectory
                ri = resample(wtilde);
                j = ri(randi(Ms));
                self.s0(:, m) = s0(:, j);
            end
        end

        %% Calculates the Backward Weights and Statistics
        % Parameters:
        %   s_p State at n+1, i.e. s[n+1] of the trajectory s[n+1:N] (Ns x
        %       1)
        %   s   Filter particles at time n, i.e. s[n] (Ns x M)
        %   w   Particle weights at time n, i.e. w[n] (1 x M)
        %   t   Time
        %   u   Control input
        %
        % Returns:
        %   W   Smoothed weights for the given backward trajectory
        %   Omega, lambda
        %       Sufficient statistcs
        %   
        % TODO: Maybe we should use the log-weights as the input.
        function [W, Omega, lambda] = calculateBackwardWeights(self, s_p, s, w, z, P, t, u)
            %% Preliminaries
            [Ns, M] = size(s);
            Nz = size(z, 1);
                        
            Is = eye(Ns);
            Iz = eye(Nz);
            Ix = eye(Ns+Nz);
            
            model = self.model;
            Omega_hat = self.Omega_hat;
            lambda_hat = self.lambda_hat;
            
            %% Calculations
            Omega = zeros(Nz, Nz, M);
            lambda = zeros(Nz, M);
            W = zeros(1, M);            
            for m = 1:M
                % 3a)
                if strcmp(self.type, 'hierarchical')                       % TODO: This selector should be modified somehow
                    % TODO: there are possibly still bugs in this branch,
                    %       related to the bugs that have been squashed 
                    %       below.
                    %% Hierarichical model
                    logZ = log(model.px_eval(s_p, s, t, u));               % TODO: Maybe the models should return the log-weight as well

                    % TODO: The calculations below don't depend on s[n] and
                    %       are thus the same for all s[n]^(i) => We only 
                    %       need to calculate these once. Should be moved
                    %       outside the loop thus.
                    % Get all the matrices / vectors from the model (in the
                    % original article's notation)
                    A = model.Al(s_p, t, u);
                    f = model.fl(s_p, t, u);
                    F = eye(Nz);   % TODo: Update as below
                    
                    % Note: These are striclty M[n+1] and m[n+1] but for
                    % convenience we call them Mt and mt anyway
                    Mt = F'*Omega_hat(:, :, m)*F + Is;
                    mt = lambda_hat(:, m) - Omega_hat(:, :, m)*f;
                    
                    % Actual calculation of the backward statistics
                    K = A'*(Iz - Omega_hat(:, :, m)*F/Mt*F');
                    Omega(:, :, m) = K*Omega_hat(:, :, m)*A;
                    lambda(:, m) = K*mt;
                elseif strcmp(self.type, 'mixed')
                    %% Mixed model
                    % Get all the matrices / vectors from the model (in the
                    % original article's notation)
                    Q = model.Qn(s(:, m), t, u);
                    g = model.fn(s(:, m), t, u);
                    B = model.An(s(:, m), t, u);
                    G = zeros(Ns, Ns + Nz);
                    G(:, model.in) = Is;
                    f = model.fl(s(:, m), t, u);
                    A = model.Al(s(:, m), t, u);
                    %F = eye(Nz);
                    F = zeros(Nz, Ns + Nz);
                    F(:, model.il) = Iz;
                    
                    % Gram-Schmidt orthogonalization
                    Qbar = Ix - G'/Q*G;
                    fbar = f + F*G'/Q*(s_p - g);
                    Abar = A - F*G'/Q*B;
                    
                    % Actual calculation of the backward statistics
                    mt = lambda_hat(:, m) - Omega_hat(:, :, m)*fbar;
                    Mt = Qbar*F'*Omega_hat(:, :, m)*F*Qbar + Ix;
                    Psi = Qbar/Mt*Qbar;
                    tau = ( ...
                        wnorm(s_p-g, Q\Is) ...
                        + fbar'*Omega_hat(:, :, m)*fbar ...
                        - 2*lambda_hat(:, m)'*fbar ...
                        - (F'*mt)'*Psi*(F'*mt) ...
                    );
                    
                    logZ = -1/2*(log(det(Q)) + log(det(Mt)) + tau);
                    K = Abar'*(Iz - Omega_hat(:, :, m)*F*Psi*F');
                    Omega(:, :, m) = K*Omega_hat(:, :, m)*Abar + B'/Q*B;
                    lambda(:, m) = K*mt + B'/Q*(s_p-g);
                else
                    error('Unknown model type.');
                end
                
                % 3b)
                zbar = z(:, m);
                Gamma = chol(P(:, :, m)).';
                Lambda = Gamma'*Omega(:, :, m)*Gamma + Iz;
                epsilon = Gamma'*(lambda(:, m) - Omega(:, :, m)*zbar);
                eta = ( ...
                    zbar'*Omega(:, :, m)*zbar ...
                    - 2*lambda(:, m)'*zbar ...
                    - epsilon'/Lambda*epsilon ...
                );
            
                % 3c)
                W(m) = log(w(m)) + logZ - 1/2*log(det(Lambda)) - 1/2*eta;
            end
            
            % 3c) Normalize the weights
            W = exp(W-max(W));
            W = W/sum(W);
        end
        
        %% Smoothing of the Linear States
        function smoothLinearStates(self, y, t, u)
            ss = self.ss;
            [~, Ms, N] = size(ss);
            Nz = size(self.z, 1);
            Iz = eye(Nz);
            model = self.model;
            zs = zeros(Nz, Ms, N);
            Ps = zeros(Nz, Nz, Ms, N);
            
            %% Process All Trajectories
            for m = 1:Ms
                s = self.s0(:, m);
                z = model.m0(model.il);
                P = model.P0(model.il, model.il);

                for n = 1:N
                    %% Hierarchical Model: Get s[n]
                    if strcmp(self.type, 'hierarchical')
                        s = self.ss(:, m, n);
                    end
                    
                    %% KF Prediction
                    % Move from n-1 to n
                    fl = model.fl(s, t, u);
                    Al = model.Al(s, t, u);
                    Ql = model.Ql(s, t, u);
                    zp = fl + Al*z;
                    Pp = Al*P*Al' + Ql;
                    
                    %% Mixed Model: Get s[n]
                    if strcmp(self.type, 'mixed')
                        s = self.ss(:, m, n);
                    end
                    
                    %% KF Measurement Update
                    h = model.h(s, t, u);
                    C = model.C(s, t, u);
                    R = model.R(s, t, u);
                    S = C*Pp*C' + R;
                    K = (Pp*C')/S;
                    z = zp + K*(y(:, n) - h - C*zp);
                    P = Pp - K*S*K';
                    
                    %% Smoothed Means & Covariances 
                    Ps(:, :, m, n) = (P\Iz + self.Omega(:, :, m, n))\Iz;
                    zs(:, m, n) = Ps(:, :, m, n)*(P\z + self.lambda(:, m, n));
                end                
            end
            
            %% Point Estimate
            % Simply the average over all particles
            zhat_s = squeeze(mean(zs, 2));
            P_s = squeeze(mean(Ps, 3));

            %% Store
            self.zs = zs;
            self.Ps = Ps;
            self.zhat_s = zhat_s;
            self.P_s = P_s;
        end
    end
end
