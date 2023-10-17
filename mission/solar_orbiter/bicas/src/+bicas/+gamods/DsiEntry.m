%
% Class that represents the GA "MODS" information for one DSI.
%
% MUTABLE.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef DsiEntry < handle
    % PROPOSAL: Assert than all comment strings (separate ones in VerionEntry)
    % are unique.

    
    
    %#####################
    %#####################
    % INSTANCE PROPERTIES
    %#####################
    %#####################
    properties(SetAccess=private, GetAccess=private)
        versionEntryCa
    end



    %#########################
    %#########################
    % PUBLIC INSTANCE METHODS
    %#########################
    %#########################
    methods(Access=public)



        function obj = DsiEntry()
            obj.versionEntryCa = {};
        end



        function add_version_entry(obj, ve)
            assert(isa(ve, 'bicas.gamods.VersionEntry'))

            % ASSERTIONS
            if ~isempty(obj.versionEntryCa)
                vePrev = obj.versionEntryCa{end};

                dateStrCa = cellfun(@(ca) ca.dateStr, ...
                    obj.versionEntryCa, 'UniformOutput', false);
                % bicasVersionStrCa = cellfun(@(ca) ca.bicasVersionStr, ...
                %     obj.versionEntryCa, 'UniformOutput', false);

                % Do not reuse earlier date.
                assert(~ismember(ve.dateStr, dateStrCa))

                % Date comes after previous one.
                % NOTE: issorted() does not check for equality for strings.
                assert(issorted({vePrev.dateStr, ve.dateStr}))

                % Do not reuse earlier BICAS version
                % ----------------------------------
                % IMPLEMENTATION NOTE: Would in principle want to test for not
                % reusing earlier BICAS version, but since there have been L3
                % deliveries made with the same BICAS version (despite minor
                % BICAS modifications), we can not use this check.                
                % assert(~ismember(ve.bicasVersionStr, bicasVersionStrCa))
                
                % NOTE: Can not easily assert sorting of BICAS version strings.
            end

            obj.versionEntryCa{end+1, 1} = ve;
        end



        function gaModsStrCa = get_MODS_strings_CA(obj)
            gaModsStrCa = cellfun(...
                @(ve) (ve.get_str()), obj.versionEntryCa(:), ...
                'UniformOutput', false);
        end



    end

end
