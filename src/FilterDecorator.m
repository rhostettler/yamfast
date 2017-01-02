classdef (Abstract) FilterDecorator < handle
    % Decorator interface definition
    %
    % DESCRIPTION
    %   Decorators add small functionalities to a filter, for example
    %   operations such as imposing constraints on the state or similar.
    %   This abstract class defines the general interface for decorators
    %   and provides dummy-implementations (that do not perform any
    %   operations). Actual decorators can then be derived from this class
    %   and re-implement the necessary functions (while using the default
    %   implementations of this abstract class where no operations are
    %   needed).
    % 
    % PROPERTIES
    %   None.
    %
    % METHODS
    %   timeUpdateHook(filter)
    %       Method that is called *after* a time update is performed by the
    %       filter.
    %
    %   measurementUpdateHook(filter)
    %       Method that is called *after* a measurement updated was
    %       performed by the filter.
    %
    % SEE ALSO
    %   NormConstraintDecorator
    % 
    % VERSION
    %   2016-12-22
    %
    % AUTHOR
    %   Roland Hostettler <roland.hostettler@aalto.fi>
        
    %% Methods
    methods (Access = public)
        %% Dummy Time Update Hook Implementation
        function timeUpdateHook(self, filter)
            % nop
        end
        
        %% Dummy Measurement Update Hook Implementation
        function measurementUpdateHook(self, filter)
            % nop
        end
    end
end
