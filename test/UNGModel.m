classdef UNGModel < AWGNModel
    % 
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
    %   2016-
    % 
    % AUTHORS
    %   Roland Hostettler <roland.hostettler@aalto.fi>   
    
    % Properties
    properties
    end

    % Methods
    methods (Access = public)
        function self = UNGModel(Q, R)
            %% Defaults
            if nargin == 0 || isempty(Q)
                Q = 0.071^2;
            end
            if nargin <= 1 || isempty(R)
                R = 0.1;
            end
            
            %% Dynamics
            f = @(x, t, u) 0.5*x + 25*x./(1+x.^2) + 8*cos(1.2*t);
            Fx = @(x, t, u) 0.5 + (25*(1+x.^2) - 50*x.^2)./(1+x.^2).^2;
            
            %% Measurement Model
            g = @(x, t, u) x.^2/20;
            Gx = @(x, t, u) x/10;

            %% Initial State
            self@AWGNModel(f, Fx, Q, g, Gx, R);
            self.m0 = 0;
            self.P0 = 1;
        end
    end
end
