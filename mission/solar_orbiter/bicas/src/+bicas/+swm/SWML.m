%
% Singleton class that stores (after having "built" it) an unmodifiable data
% structure that represents which and how s/w modes are CURRENTLY VISIBLE to the
% user. What that data structure contains thus depends on
% -- current pipeline: RODP, ROC-SGSE
% -- whether support for L1 input datasets is enabled or not.
%
% Data here
% -- is intended to (1) add metadata to the production functions for
%       (a) the caller interface
%       (b) potentially future help text
%       (c) the s/w descriptor
% -- contains variables to make it possible to match information for input and
%    output datasets here, with that of bicas.proc' production functions.
%
%
% IMPLEMENTATION NOTES
% ====================
% The class essentially consists of one large struct, and a constructor that
% builds it. The large data struct contains many parts which are similar but not
% the same. To do this, much of the data is "generated" with hard-coded strings
% (mostly the same in every iteration), in which specific codes/substrings are
% substituted algorithmically (different in different iterations). To avoid
% mistakes, the code uses a lot of assertions to protect against mistakes, e.g.
% -- algorithmic bugs
% -- mistyped hard-coded info
% -- mistakenly confused arguments with each other.
% Assertions are located at the place where "values are placed in their final
% location".
% 
% NOTE: To implement compatibility with L1 input datasets, the code must be able
% to handle
% -- changing input dataset levels: L1 (unofficial support), L1R (official
%    support). It implements support for L1 input datasets via separate S/W
%    modes.
% 
%
% RATIONALE
% =========
% -- Should decrease the amount of overlapping hard-coded information to e.g.
%    reduce risk of mistakes, reduce manual work when verifying updates.
% -- Having one big, albeit somewhat redundant data structure should make the
%    interface to the rest of BICAS relatively future-proof, in the face of
%    future updates.
% -- Useful for expected future bias current datasets.
% -- Possible need for backward compatibility.
%
%
% NOTE
% ====
% There is no s/w mode for generating VHT datasets. They are therefore not
% represented here.
%
%
% DEFINITIONS
% ===========
% SIP = "Specific Input Parameters" (RCS ICD).
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2019-07-31
%
classdef SWML
    % PROPOSAL: Pick SWD name/descriptions from master CDFs.
    % PROPOSAL: Obtain output dataset level from production function metadata?!!
    %
    % PROPOSAL: Always produce all possible s/w modes (both pipelines, incl. L1),
    %           then filter out the undesired ones using internal metadata for
    %           every S/W mode.
    %
    % PROPOSAL: Use PF = prodFunc, production function
    % PROPOSAL: Same input CDF can have multiple DATASET_IDs, but only one is
    %           shown in the s/w descriptor.
    %   PRO: Can handle old datasets with ROG-SGSE DATASET_IDs, and otherwise
    %        only use RODP DATASET_IDs.
    %
    % TODO-DEC: Which arguments should SWML production functions (function handles in
    %           an instance of bicas.swm.SWML) have?
    %   NOTE: The arguments needed by the underlying production functions
    %         varies, but the arguments returned by bicas.swm.SWML must be the same.
    %   NOTE: produce_L1R_to_L2_LFR/TDS() are used for multiple s/w modes with some
    %         arguments hard-coded differently for different s/w modes (input & output DATASET_IDs).
    %   NOTE: SWM/underlying production functions can receive argument values via
    %       (1) bicas.swm.SWML (constructor), or (2) the call in execute_SWM.
    %   PROPOSAL: All arguments which are known at the time bicas.swm.SWML
    %       constructor is called, should receive values there.
    %       ==> ~As many as possible.
    %       CON: bicas.swm.SWML not really meant to set production function arguments.
    %       CON: Makes bicas.swm.SWML harder to initialize (outside of BICAS).
    %   PROPOSAL: All arguments which are different for different (underlying) production
    %             functions. ==> As few as possible.
    %   Ex: SETTINGS, L, rctDir, NsoTable
    %
    % PROPOSAL: Abolish this class. Only need SWM class + object array/list.
    % PROPOSAL: Better class name
    %   PROPOSAL: SwmSet
    %       PRO: There is no inherent ordering of SWMs.
    %   PROPOSAL: SwmList
    %       CON: List implies ordering.
    %       CON: Conflicts with variable naming for plain list of SWM objects.



    % PUBLIC, IMMUTABLE
    properties(SetAccess=immutable)
        
        % NOTE: Implicit that it is a list of s/w modes (not in name). Note that
        % it is a public property.
        List
    end
    
    
    
    methods(Access=public)
        
        
        
        % Constructor
        %
        % ARGUMENTS
        % =========
        %
        % IMPLEMENTATION NOTE: The constructor used to be written so that it was
        % easy to disable S/W modes with L2R input datasets (for backward
        % compatibility). That functionality has now been now removed, although
        % the implementation has not been entirely updated to take advantage of
        % this (not simplified of this).
        % 
        function obj = SWML(SwmList)
            assert(isvector(SwmList))
            assert(isa(SwmList, 'bicas.swm.SWM'))
            irf.assert.castring_set({SwmList(:).cliOption})
            
            obj.List = SwmList;
        end



        function swm = get_SWM(obj, swmCliOption)
            i = find(strcmp(swmCliOption, {obj.List(:).cliOption}));
            irf.assert.scalar(i)
            swm = obj.List(i);
        end
        
        
        
    end    % methods(Access=public)


    
end
