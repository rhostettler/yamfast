classdef MixedRBPF < handle
    % Rao-Blackwellized particle filter for mixed linear/non-linear models
    % 
    % DESCRIPTION
    %   Rao-Blackwellized particle filter for mixed linear/non-linear
    %   models according to [1].
    %
    %   
    %
    % PROPERTIES
    % 
    %
    % METHODS
    %
    %
    % REFERENCES
    %   [1] T. Schön, F. Gustafsson, P.-J. Nordlund, "Marginalized Particle
    %       Filters for Mixed Linear/Nonlinear State-Space Models", IEEE
    %       Transactions on Signal Processing", vol. 53, no. 7, July 2005
    %
    % SEE ALSO
    %   RBFFBSi
    %
    % VERSION
    %   2017-01-17
    % 
    % AUTHORS
    %   Roland Hostettler <roland.hostettler@aalto.fi>
    
     %#ok<*PROPLC>
    
    %% Properties
    properties (Access = public)
        % The model
        model;
        
        % Non-linear states & weights
        s;
        w;
        
        % Linear states & covariance
        z;
        P;
        
        % Resampling threshold
        Mt;
    end
    
    %% Pseudo-Properties
    properties (Dependent)
    end

    %% Methods
    methods (Access = public)
        %% Initialization
        function self = MixedRBPF(model, M)
            self.model = model;
            m0 = model.m0;
            P0 = model.P0;
            Ns = length(model.in);
            self.s = m0(model.in) + chol(P0(model.in, model.in)).'*randn(Ns, M);
            self.w = 1/M*ones(1, M);
            self.z = m0(model.il)*ones(1, M);
            self.P = repmat(P0(model.il, model.il), [1, 1, M]);
            self.Mt = M/3;
        end
        
        %% Filter Update
        function xhat_f = update(self, y, t, u)
            model = self.model;
            s = self.s;
            w = self.w;
            z = self.z;
            P = self.P;

            [Ns, M] = size(s);
            Nz = size(z, 1);
            Is = eye(Ns);
            Iz = eye(Nz);
            
            z_p = zeros(size(z));
            P_p = zeros(size(P));
            
            for m = 1:M
                %% Dynamics from n-1 to n
                % Get matrices
                fn = model.fn(s(:, m), t, u);
                An = model.An(s(:, m), t, u);
                Gn = Is;
                Qn = model.Qn(s(:, m), t, u);
                
                fl = model.fl(s(:, m), t, u);
                Al = model.Al(s(:, m), t, u);
                Gl = Iz;
                Ql = model.Ql(s(:, m), t, u);
                Qnl = model.Qnl(s(:, m), t, u);

                % Orthogonalization
                Albar = Al - Gl*Qnl'/(Gn*Qn)*An;
                Qlbar = Ql - Qnl'/Qn*Qnl;

                %% Draw new Particles
                mu = fn + An*z(:, m);
                Sigma = An*P(:, :, m)*An' + Gn*Qn*Gn';
                s(:, m) = mu + chol(Sigma).'*randn(Ns, 1);
            
                %% KF Prediction
                Nt = An*P(:, :, m)*An' + Gn*Qn*Gn';
                Lt = Albar*P(:, :, m)*An'/Nt;
                v = s(:, m) - fn;
                z_p(:, m) = ( ...
                    Albar*z(:, m) + Gl*Qnl'/(Gn*Qn)*v ...
                    + fl + Lt*(v - An*z(:, m)) ...
                );
                P_p(:, :, m) = ( ...
                    Albar*P(:, :, m)*Albar' + Gl*Qlbar*Gl' - Lt*Nt*Lt' ...
                );
                
                %% KF Measurement update
                C = model.C(s(:, m), t, u);
                R = model.R(s(:, m), t, u);
                h = model.h(s(:, m), t, u);
                Mt = C*P_p(:, :, m)*C' + R;
                Kt = P_p(:, :, m)*C'/Mt;
                z(:, m) = z_p(:, m) + Kt*(y - h - C*z_p(:, m));
                P(:, :, m) = P_p(:, :, m) - Kt*Mt*Kt';
                
                %% Weights
                % TODO: Implement as log instead
                mu = h + C*z_p(:, m);
                Sigma = C*P_p(:, :, m)*C' + R;
                % w(:, m) = log(w(:, m)) + logmvnpdf(y, mu, Sigma);
                w(:, m) = w(:, m)*mvnpdf(y, mu, Sigma);
            end
            
            %% Point Estimate and Resampling
            % Normalize the weights
            %w = exp(w-max(w));
            w = w/sum(w);
            
            % MMSE
            xhat_f = zeros(Ns + Nz, 1);
            xhat_f(model.in) = s*w';
            xhat_f(model.il) = z*w';
            
            % Resampling
            ess = 1/sum(w.^2);
            if ess <= self.Mt
                ri = resample(w);
                s = s(:, ri);
                z = z(:, ri);
                P = P(:, :, ri);
                w = 1/M*ones(1, M);
            end
            
            %% Store
            self.s = s;
            self.w = w;
            self.z = z;
            self.P = P;
        end
    end
end
