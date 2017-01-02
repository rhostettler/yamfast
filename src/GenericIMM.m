classdef GenericIMM < handle
    % Interacting multiple model
    % 
    % DESCRIPTION
    %   A generic ("container") model for interacting multiple models. This 
    %   model class encapsules other models and the corresponding switching
    %   and initial probabilities, as well as the state mapping between the
    %   models.
    %
    % PROPERTIES
    %   models (r, cell array, 1xM)
    %       The actual storage container of the models.
    %
    %   mappings (r, cell array, 1xM)
    %       Stores the state mapping of the individual models' state vector
    %       to the global, common state vector.
    %
    %   p_ij (r/w, matrix, MxM)
    %       Matrix of transition probabilities from model i to model j.
    %       Rows and columns must sum to one.
    %
    %   mu_i0 (r/w, vector, Mx1)
    %       Initial probabilities of the different models.
    %
    %   M (r)
    %       No. of models
    %
    %   Nx (r)
    %       Dimension of the common state vector.
    %
    % METHODS
    %   addModel(model, mapping)
    %       Add a new model with the specified state mapping to the IMM.
    %
    % SEE ALSO
    %   GenericIMMFilter
    %
    % VERSION
    %   2017-01-02
    % 
    % AUTHORS
    %   Roland Hostettler <roland.hostettler@aalto.fi>   
    
    %% Properties
    % Read/write
    properties (Access = public)
        % Initial probabilities
        mu_i0 = [];
        
        % Transition probabilities
        p_ij = [];
    end
    
    % Read-only
    properties (GetAccess = public)
        % Models
        models = {};
        
        % State mappings
        mappings = {};
        
        % Total no. of models
        M = 0;
        
        % Dimension of common state
        Nx = 0;
    end

    %% Public Methods
    methods (Access = public)
        function addModel(self, model, mapping)
            self.M = self.M+1;
            self.models{self.M} = model;
            self.mappings{self.M} = mapping;            
            self.Nx = max([self.Nx; mapping(:)]);
        end
        
        function [x_p, F] = f(self, x, q, t, u, i)
            [x_p, F] = self.models{i}.f(x, q, t, u);
        end
        
        function [y, G] = g(self, x, r, t, u, i)
            [y, G] = self.models{i}.g(x, r, t, u);
        end
        
        function setTime(self, t)
            for m = 1:self.M
                self.models{m}.t = t;
            end
        end
    end
end
