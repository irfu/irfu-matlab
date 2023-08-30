%
% Class whose instances represent one loaded dataset (CDF file).
%
% NOTE: Not to be confused with bicas.swm.InputDataset.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-09-09
%
classdef InputDataset
    % PROPOSAL: Make all properties/instance variables immutable. Modify by
    %           creating modified copy of class instance.



    %#####################
    %#####################
    % INSTANCE PROPERTIES
    %#####################
    %#####################
    
    %====================
    % Mutable properties
    %====================
    properties
        % Struct with zVariable values.
        % IMPLEMENTATION NOTE: Not immutable since some values are normalized
        % during processing (2x6x L1/L1R-->L2).
        Zv
    end
    
    %======================
    % Immutable properties
    %======================
    properties(SetAccess=immutable)
    
        % Struct with zVariable fill values.
        ZvFv
        
        % Struct with global attributes.
        Ga
        
        filePath
        
    end    % properties(SetAccess=immutable)



    %#########################
    %#########################
    % PUBLIC INSTANCE METHODS
    %#########################
    %#########################
    methods(Access=public)
        
        
        
        function obj = InputDataset(Zv, ZvFv, Ga, filePath)
            obj.Zv       = Zv;
            obj.ZvFv     = ZvFv;
            obj.Ga       = Ga;
            obj.filePath = filePath;
        end
        
        
        
    end    % methods(Access=public)



end
