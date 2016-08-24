classdef (Abstract) GenericStateSpaceModel < handle
    % GenericStateSpaceModel Interface for generic state space models
    % 
    % DESCRIPTION
    %   This (abstract) class defines the interface for generic state-space
    %   models of the form
    %
    %       x(k) = f(x(k-1), q(k-1), t(k))
    %       y(k) = g(x(k), r(k), t(k))
    %
    %   where q(k) and r(k) are the process- and measurement noises with
    %   covariance Q and R, respectively, f(.) is the state transition
    %   function, and g(.) is the observation function. t(k) is the time at
    %   sampling instant k.
    %
    %   Classes implementing this model type need to implement the methods 
    %   'f' and 'g' according to the above definition.
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
    %   2016-03-01
    % 
    % AUTHORS
    %   Roland Hostettler <roland.hostettler@aalto.fi>  
    
    % TODO
    %   * This is actually not exactly the most generic case (especially
    %     w.r.t. the noise variables and initialization). I should try to
    %     generalize it even more.
    %   * Also, Q and R could be time dependent.
    %   * All functions/values could be parameter dependent.
    
    % Properties
    properties
        % Initial state and covariance
        m0;
        P0;
        
        % Process noise covariance
        Q;
        
        % Measurement noise covariance
        R;
    end

    % Abstract methods, these need to be implemented by derived classes
    methods (Abstract)
        %% State Dynamics
        f(self, x, u, v, t);
        
        %% Measurement Function
        g(self, x, u, w, t);
    end
end
