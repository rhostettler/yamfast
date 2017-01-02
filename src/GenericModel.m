classdef (Abstract) GenericModel < handle
    % Model interface for generic state space models
    % 
    % DESCRIPTION
    %   This (abstract) class defines the interface for generic state-space
    %   models of the form
    %
    %       x[n] = f(x[n-1], q[n], t[n], u[n])
    %       y[n] = g(x[n], r[n], t[n], u[n])
    %
    %   where q[n] and r[n] are the process- and measurement noises, 
    %   respectively, f(.) is the state transition function, and g(.) is 
    %   the observation function. t[n] is the time at sampling instant n
    %   and u[n] is a possible deterministic control input.
    %
    %   This class defines the most basic interface for such models.
    %   Classes implementing this model type need to implement the
    %   methods listed in the section 'ABSTRACT METHODS' below.
    %
    % PROPERTIES
    %   t (r/w)
    %       The time of the last update.
    %
    % ABSTRACT METHODS
    %   x0 = px0_rand(M)
    %       Random state generator for the initial state where M is the
    %       number of states to generate.
    %
    %   [xp, Fx, Fq] = f(x, q, t, u)
    %       State transition function where xp is the predicted state, and 
    %       Fx and Fq the Jacobians of f(.) with respect to x and q,
    %       respectively.
    %
    %   xp = px_rand(x, t, u)
    %       Function to generate random states from the state transition
    %       density p(x[k] | x[k-1]).
    %
    %   px = px_eval(xp, x, t, u)
    %       Function to evaluate the state transition density.
    %
    %   [yp, Gx, Gr] = g(x, r, t, u)
    %       Measurement function where yp is the predicted measurement and
    %       Gx and Gr are the Jacobians of g(.) with respect to x and r,
    %       respectively.
    %
    %   py = py_eval(x, t, u)
    %       Measurement likelihood.
    %
    % SEE ALSO
    %   AWGNModel, LGSSModel, MixedCLGSSModel
    %
    % VERSION
    %   2017-01-02
    % 
    % AUTHORS
    %   Roland Hostettler <roland.hostettler@aalto.fi>  
        
    %% Properties
    properties
        % Time of the last update for time-dependent models. This needs to 
        % be stored by the main application after each simulation/filter
        % update.
        t = 0;
    end

    %% Abstract Methods
    % These need to be implemented by child-classes
    methods (Abstract)
        px0_rand(self, M);            % Initial state generator
        f(self, x, q, t, u);          % State Dynamics
        px_rand(self, x, t, u);       % State transition generator
        px_eval(self, xp, x, t, u);   % State transition evaluation
        g(self, x, r, t, u);          % Measurement function
        py_eval(self, y, x, t, u);    % Measurement likelihood
    end
end
