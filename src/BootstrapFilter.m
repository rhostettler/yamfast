classdef BootstrapFilter < handle
    % Bootstrap Particle Filter
    % 
    % DESCRIPTION
    %   A standard SIR bootstrap particle filter targeting the joint
    %   filtering distribution p(x_{1:t} | y_{1:t}) for state-space models 
    %   of the type
    %
    %       x[n] ~ p(x[n] | x[n-1])
    %       y[n] ~ p(y[n] | x[n])
    %
    %   The filter uses effective sample size-based resampling by default
    %   with a default threshold of M/3 as well as systematic resampling.
    %
    % PROPERTIES
    %   model
    %       The model the filter operates on. See GenericModel for details
    %       on the model requirements.
    %
    %   x, w
    %       The particles and their weights at time t.
    %
    %   Mt  Resampling threshold.
    %
    %   xhat, xmap, P
    %       Posterior mean, mode, and covariance at time t.
    %
    % METHODS
    %   BootstrapFilter(model, M)
    %       Initialize the bootstrap filter with M particles for the given
    %       model.
    %
    %   [xhat, P] = update(y, t, u)
    %       Do a complete update with measurement y, time t, and control
    %       input u.
    %
    % SEE ALSO
    %   GenericModel
    %
    % VERSION
    %   2017-01-02
    % 
    % AUTHORS
    %   Roland Hostettler <roland.hostettler@aalto.fi>   
    
    % TODO
    %   * Should we derive this from a common filter class?
    %   * Move resampling out of the class
    %   * Add bells and whistles like jittering and stuff
    %   * Add input check function
    %   * Add possibility for fast propagation and evaluation
    %   * Add diagnostics (particle degeneracy, etc)
    %   * Store (and resample) whole particle trajectories
    
    % Supress mlint warnings
    %#ok<*PROPLC>
    
    %% Properties
    properties (Access = public)
        % Model
        model;
        
        % Particles & their weights
        x;
        w;
        
        % Resampling threshold, defaults to M/3
        Mt;
        
        % Posterior mean and mode and covariance
        xhat;
        xmap;
        P;
    end

    %% Methods
    methods (Access = public)
        %% Initialization
        function self = BootstrapFilter(model, M)
            self.model = model;
            self.Mt = M/3;
            self.x = model.px0_rand(M);
            self.w = 1/M*ones(1, M);
        end
        
        %% Filter Update
        function [xhat, P]  = update(self, y, t, u)
            model = self.model;
            w = self.w;
            x = self.x;
            M = size(x, 2);
            
            % Propagate
            for m = 1:M
                x(:, m) = model.px_rand(x(:, m), t, u);
                w(:, m) = w(:, m)*model.py_eval(y, x(:, m), t, u);
            end
            w = w/sum(w);
            
            % Estimate
            xhat = x*w';
            P = 0;
            for m = 1:M
                P = P + w(m)*((x(:, m) - xhat)*(x(:, m) - xhat)');
            end
            [~, imap] = max(w);
            self.xmap = x(:, imap);
            
            % Resample
            ess = 1/sum(w.^2);
            if ess < self.Mt
                ir = resample(w);
                x = x(:, ir);
                w = 1/M*ones(1, M);
            end
                        
            % Store
            self.x = x;
            self.w = w;
            self.xhat = xhat;
            self.C = P;
        end
    end
    
    %% Helper Methods
    methods (Access = private)
        %% Systematic Resampling
        function ri = resample(self, w)
            M = length(w);
            ri = zeros(1, M);
            i = 0;
            u = 1/M*rand();
            for j = 1:M
                Ns = floor(M*(w(j)-u)) + 1;
                if Ns > 0
                    ri(i+1:i+Ns) = j;
                    i = i + Ns;
                end
                u = u + Ns/M - w(j);
            end
            ri = ri(randperm(M));
        end
    end
end
