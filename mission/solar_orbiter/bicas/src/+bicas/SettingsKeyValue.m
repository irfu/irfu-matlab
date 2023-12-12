%
% Class which encapsulates the information stored for one settings key in
% bicas.Settings.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef SettingsKeyValue
% PROPOSAL: Better name
%   ~Settings
%       CON: Want to avoid if moving class to package for settings.
%   ~Entry
%   ~Key, Value, Pair
%   Key
%       PRO: There is one key per object but multiple values.
%   SettingsEntry
%       CON: Two-letter abbreviation.
%   SettingsKeyEntry
%   SettingsKeyValueEntry = SKVE
%       NOTE: Can drop "Settings" if moving class to settings package.



    %#####################
    %#####################
    % INSTANCE PROPERTIES
    %#####################
    %#####################
    properties(SetAccess=private, GetAccess=public)
        % Column cell array of values.
        %
        % NOTE: Must be CA of values since the values might (1) be strings, or
        % (2) arrays of varying sizes.
        valuesCa

        % Column cell array.
        valueSourcesCa
    end



    %#########################
    %#########################
    % PUBLIC INSTANCE METHODS
    %#########################
    %#########################
    methods(Access=public)



        function obj = SettingsKeyValue(initialValue, valueSource)
            assert(ischar(valueSource))

            obj.valuesCa       = {initialValue};
            obj.valueSourcesCa = {valueSource};
        end



        function obj = override(obj, value, valueSource)
            assert(ischar(valueSource))
            assert(~ismember(valueSource, obj.valueSourcesCa))

            % NOTE: Creating column arrays.
            obj.valuesCa{end+1, 1}       = value;
            obj.valueSourcesCa{end+1, 1} = valueSource;
        end
        
        
        
        function n = N_values(obj)
            n = numel(obj.valuesCa);
        end



    end    % methods(Access=public)



end
