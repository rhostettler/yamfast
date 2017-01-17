classdef GenericIMMFilter < handle
    % Generic interacting multiple model filter
    % 
    % DESCRIPTION
    %   (Gaussian) Filter for interacting multiple models according to [1]. 
    %   Does all the mixing and pruning and can use any Gaussian filter.
    %
    % PROPERTIES
    %   m, P, m_p, P_p
    %       Combined posterior mean and covariance, as well as predicted
    %       mean and covariance.
    %
    %   mu_i
    %       Model probabilities.
    %
    %   m_i, P_i, m_pi, P_pi
    %       Individual posterior means, covariances, as well as individual
    %       predicted means and covariances.
    %
    % METHODS
    %   [m, P] = update(y, t, u)
    %       Filter update.
    %
    % SEE ALSO
    %   GenericIMM
    %
    % VERSION
    %   2017-01-02
    % 
    % AUTHORS
    %   Roland Hostettler <roland.hostettler@aalto.fi>   
    
    % TODO: Should move the predictive likelihood to the underlying filter
    %       instead of here, then we can use it with other types of assumed
    %       density filters too.
    
    % Suppress mlint warnings
    %#ok<*PROPLC>
    
    %% Properties
    properties (Access = public)        
        % Combined mean, covariance and predicted mean and covariance
        m;
        P;
        m_p;
        P_p;
        
        % Model probabilities
        mu_i;

        % Predicted state and covariance for each model (1xM cell array)
        m_i;
        P_i;
        m_pi;
        P_pi;        
    end
    
    %% Private Properties
    properties (Access = private)
        % Stores the filter
        filter;
        
        % Stores the IMM
        imm;
    end

    %% Methods
    methods (Access = public)
        function self = GenericIMMFilter(model, filter)
            self.filter = filter;
            self.imm = model;
            M = model.M;
            
            % Initialize the storage structures
            self.m_i = cell(1, M);
            self.P_i = cell(1, M);
            self.mu_i = zeros(M, 1);
            self.m_pi = cell(1, M);
            self.P_pi = cell(1, M);
            
            % Initialize with the intial values
            for i = 1:M
                self.m_i{i} = model.models{i}.m0;
                self.P_i{i} = model.models{i}.P0;
            end
            
            % Get the initial probabilities
            self.mu_i = model.mu_i0;
        end
                
        function [m, P] = update(self, y, t, u)
            %% Commonly Used Variables
            imm = self.imm;   % The model
            M = imm.M;        % M: No. of models
            Nx = imm.Nx;      % Nx: Largest state dimension
            
            % 'mappings' -> index vector with indices of the state
            % components of the i-th model in the combined state vector
            mappings = imm.mappings;
            filter = self.filter;
            
            %% Interaction
            % Mixing probabilites
            mu_ij = zeros([M, M]);
            c_j = zeros(M, 1);
            for j = 1:M
                c_j(j) = imm.p_ij(:, j).'*self.mu_i;
                mu_ij(:, j) = 1/c_j(j)*imm.p_ij(:, j).*self.mu_i;
            end
            
            % Mixed input mean
            m_0j = cell(1, M);
            for j = 1:M
                m_0j{j} = zeros(Nx, 1);
                for i = 1:M
                    m_0j{j}(mappings{i}) = m_0j{j}(mappings{i}) + mu_ij(i, j)*self.m_i{i};
                end
            end

            % Mixed input covariance
            P_0j = cell(1, M);
            for j = 1:M
                P_0j{j} = zeros(Nx, Nx);
                for i = 1:M
                    P_0j{j}(mappings{i}, mappings{i}) = ( ...
                        P_0j{j}(mappings{i}, mappings{i}) ...
                        + mu_ij(i, j)*( ...
                            self.P_i{i} ...
                            + (self.m_i{i}-m_0j{j}(mappings{i}))*(self.m_i{i}-m_0j{j}(mappings{i}))' ...
                        ) ...
                    );
                end
            end
            
            %% Filter
            lambda_j = zeros(M, 1);
            for j = 1:M
                % Set the currently active model, mean, and covariance
                filter.model = imm.models{j};
                filter.m = m_0j{j}(mappings{j});
                filter.P = P_0j{j}(mappings{j}, mappings{j});
                
                % Update (using the standard UKF's update function)
                filter.update(y, t, u);
                
                % Calculate the likelihood through the innovation
                lambda_j(j) = mvnpdf((y-filter.y_p).', zeros(1, size(y, 1)), filter.S).';
                
                % Store mean and covariance (both predicted and corrected)
                self.m_pi{j} = filter.m_p;
                self.P_pi{j} = filter.P_p;
                self.m_i{j} = filter.m;
                self.P_i{j} = filter.P;
            end
            
            % Calculate the model probabilities
            c = lambda_j'*c_j;
            self.mu_i = lambda_j.*c_j/c;
                            
            %% Combination
            % Mean
            m = zeros(Nx, 1);
            for i = 1:M
                m(mappings{i}) = m(mappings{i}) + self.mu_i(i)*self.m_i{i};
            end
            
            % Covariance
            P = zeros(Nx, Nx);
            for i = 1:M
                P(mappings{i}, mappings{i}) = ( ...
                    P(mappings{i}, mappings{i}) ...
                    + self.mu_i(i)*( ...
                        self.P_i{i} ...
                        + (self.m_i{i} - m(mappings{i}))*(self.m_i{i} - m(mappings{i}))'...
                    ) ...
                );
            end
            
            self.m = m;
            self.P = P;
        end
    end
end
