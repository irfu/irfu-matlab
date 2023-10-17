%
% Class that represents the GA "MODS" information for one DSI.
% Contains list of instances of bicas.gamods.VersionEntry.
% 
%
% MUTABLE. HANDLE CLASS.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef DsiEntry < handle
% PROPOSAL: Assert than all comment strings (separate ones in VersionEntry)
%           are unique.
% PROPOSAL: Rename add_version_entry() --> add_GMVE().

    
    
    %#####################
    %#####################
    % INSTANCE PROPERTIES
    %#####################
    %#####################
    properties(SetAccess=private, GetAccess=private)
        gmveCa
    end



    %#########################
    %#########################
    % PUBLIC INSTANCE METHODS
    %#########################
    %#########################
    methods(Access=public)



        function obj = DsiEntry()
            obj.gmveCa = {};
        end



        function add_version_entry(obj, Gmve)
            assert(isa(Gmve, 'bicas.gamods.VersionEntry'))

            % ASSERTIONS
            if ~isempty(obj.gmveCa)
                GmvePrev = obj.gmveCa{end};

                dateStrCa = cellfun(@(ca) ca.dateStr, ...
                    obj.gmveCa, 'UniformOutput', false);
                % bicasVersionStrCa = cellfun(@(ca) ca.bicasVersionStr, ...
                %     obj.gmveCa, 'UniformOutput', false);

                % Do not reuse earlier date.
                assert(~ismember(Gmve.dateStr, dateStrCa))

                % Date comes after previous one.
                % NOTE: issorted() does not check for equality for strings.
                assert(issorted({GmvePrev.dateStr, Gmve.dateStr}))

                % Do not reuse earlier BICAS version
                % ----------------------------------
                % IMPLEMENTATION NOTE: Would in principle want to test for not
                % reusing earlier BICAS version, but since there have been L3
                % deliveries made with the same BICAS version (despite minor
                % BICAS modifications), we can not use this check.                
                % assert(~ismember(Gmve.bicasVersionStr, bicasVersionStrCa))
                
                % NOTE: Can not easily assert sorting of BICAS version strings.
            end

            obj.gmveCa{end+1, 1} = Gmve;
        end



        % Return cell array of strings to be used as value GA MODS for the
        % specified DSI.
        function gaModsStrCa = get_MODS_strings_CA(obj)
            gaModsStrCa = cellfun(...
                @(Gmve) (Gmve.get_str()), obj.gmveCa(:), ...
                'UniformOutput', false);
        end



    end



end
