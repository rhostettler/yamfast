classdef (Abstract) GenericModel < handle
    % Model Interface for generic state space models
    % 
    % DESCRIPTION
    %   This (abstract) class defines the interface for generic state-space
    %   models of the form
    %
    %       x[k] = f(x[k-1], q[k], t[k], u[k])
    %       y[k] = g(x[k], r[k], t[k], u[k])
    %
    %   where q[k] and r[k] are the process- and measurement noises, 
    %   respectively, f(.) is the state transition function, and g(.) is 
    %   the observation function. t[k] is the time at sampling instant k
    %   and u[k] is a possible deterministic control input.
    %
    %   Classes implementing this model type need to implement the
    %   following methods listed in the section 'ABSTRACT METHODS'.
    %
    %   []
    %
    %   'f' and 'g' according to the above definition.
    %
    % PROPERTIES
    %   t   The time of the last update.
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
    %   [yp, Gx, Gr] = g(x, r, t, u)
    %       Measurement function where yp is the predicted measurement and
    %       Gx and Gr are the Jacobians of g(.) with respect to x and r,
    %       respectively.
    %
    %   py = py_eval(x, t, u)
    %       Measurement likelihood.
    %
    % SEE ALSO
    %   LGSSModel, MixedLNGModel
    %
    % VERSION
    %   2016-10-19
    % 
    % AUTHORS
    %   Roland Hostettler <roland.hostettler@aalto.fi>  
        
    %% Properties
    properties
        t;                      % Model time
    end

    %% Abstract Methods
    % These need to be implemented by child-classes
    methods (Abstract)
        px0_rand(M);            % Initial state generator
        f(self, x, q, t, u);    % State Dynamics
        px_rand(x, t, u);       % State transition generator
        g(self, x, r, t, u);    % Measurement function
        py_eval(x, t, u);       % Measurement likelihood
    end
end
