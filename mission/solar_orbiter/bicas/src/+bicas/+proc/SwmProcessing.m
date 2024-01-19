%
% Abstract class which (subclass) instances represent the processing associated
% with one SWM.
%
% NOTE: Not to be confused with class bicas.swm.SoftwareMode which encapsulates
% an instance of (subclasses of) this class, plus metadata needed for the SWD.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef(Abstract) SwmProcessing
    % PROPOSAL: Naming convention for subclasses: *Swmp
    %   PRO: Shortens names of implementations.
    %       Ex:
    %         bicas.proc.L1L2.LfrSwmProcessing
    %         bicas.proc.L1L2.TdsSwmProcessing
    %         bicas.proc.L2L2.LfrDsrSwmProcessing
    %         bicas.proc.L2L3.L3OsrDsrSwmProcessing
    %   CON: Bad for TestSwmProcessing, SwmProcessing
    %       CON: Does not have to apply convention to every class.



    %###########################
    %###########################
    % ABSTRACT INSTANCE METHODS
    %###########################
    %###########################
    methods(Abstract)



        % ARGUMENTS
        % =========
        % InputsMap
        %       containers.Map with
        %       <keys>   : String defining a name of an input ("prodFuncInputKey" in
        %                  bicas.swm.SoftwareModeList).
        %       <values> : A struct with data corresponding to a CDF file
        %                  (zVariables+global attributes).
        % OutputsMap
        %       containers.Map with
        %       <keys>   : String defining a name of an output ("prodFuncOutputKey" in
        %                  bicas.swm.SoftwareModeList).
        %       <values> : A struct with data corresponding to a CDF file (zVariables).
        % --
        % Production functions should not assume/specify any particular input
        % dataset version, but should read it out from global attributes (part
        % of the PDV) if necessary.
        OutputDatasetsMap = production_function(obj, ...
            InputDatasetsMap, rctDir, NsoTable, Bso, L)
        % PROPOSAL: Rename
        %   process
        %   process_data
        %   process_method
        %
        % PROPOSAL: Remove argument rctDir, NsoTable.
        %   PRO: Not used for L2-L2, L2-L3.
        %   CON: rctDir, NsoTable are not available when constructor is called.
        %       PRO: Do not want to load NsoTable in case the SWM that actually
        %            runs does not use it.
        %           CON: The current implementation always loads NsoTable.
        %   PROPOSAL: Submit path to NSO table in (subclass) constructor.



    end    % methods(Access=public)



end
