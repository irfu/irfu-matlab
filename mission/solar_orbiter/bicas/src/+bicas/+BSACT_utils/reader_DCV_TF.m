%
% Class that collects calibration data and the corresponding metadata from BSACT. An object can and should contain
% either DCV or TF data. Metadata are collected by parsing the testlogbook*.txt files.
%
%
%
% NOTES
% =====
% NOTE: The implementation does not force one to use either DCV or TF test in a given object, one can mix them, but
% that form uf usage makes little sense and is not intended.
% NOTE: There is a slight difference between BSACT TFs (without inverted diffs).
%   2016 June: Has TF phase shift ~180 degrees at low frequencies.
%   2016 July: Has TF phase shift   ~0 degrees at low frequencies.
% NOTE: The class name is chosen to reflect the types of calibration data that it may contain in anticipation of 
% eventually creating another analogous class for other calibration data (bias current calibration data).
% 
%
%
% RATIONALE
% =========
% Loading all calibration files rather than just their fits, or the relevant subsets could be useful:
% (1) Can be used for automatically deriving smaller sets of calibration data to actually be used by BICAS.
%     This includes deriving more "refined" (more granular), smaller sets of calibration data to actually be used by
%     BICAS:
%       Example 1: Select AC (spectrum) transfer functions (TF) based on combination of "channel" and high gain/low
%       gain, rather than just one for AC low gain and AC high gain as expected.
%       Example 2: Interpolate table instead of linear function (offset & slope; y=k*x+m).
%       Example 3: Derive fits (DC) from calibration tables directly rather than using separately tabulated offsets and
%       slopes.
% (2) Can be used for automatically converting calibration data to to other file formats (CDF, XML etc).
% (3) Manually analyzing the calibration data itself (plotting; making one's own fits using other table columns;
% comparing channels, checking temperature dependence).
%
%
%
% VARIABLE NAMING CONVENTIONS / DEFINITIONS OF TERMS
% ==================================================
% BSACT             BIAS StandAlone Calibration Tables. The raw files generated at the BIAS standalone calibration, in
%                   particular taken June 2016 and July 2016.
%
% CTable:           Calibration table. Equal to the contents of one BIAS standalone calibration file.
%
% metadata          Refers to the metadata for a specific calibration table/CTable, e.g. MUX mode, input channels,
%                   AC low/high gain, diff/single, stimuli.
%
% Input channel:    One or two integers (row vector; values 1-3), signifying a single-probe "channel" or
% (inputChNbrs)     a diff-probe "channel", i.e. which signal is coming from the ANTENNA(S).
%                       Example: [2] refers to V2_LF.
%                       Example: [2,3] refers to (V2_LF-V3_LF).
%                   NOTE: Does NOT refer to the signal itself (voltage), only which channel out of several.
%                   NOTE: The probe order is always increasing.
%
% Output channel:   Integer (1-5) representing one of the "channels" at the BIAS-to-LFR/TDS boundary.
% (outputChNbr)     Also known as LFR_1...LFR_5, TDS_1...TDS_3, and BIAS_1...BIAS_5.
%                   NOTE: Does NOT refer to the signal itself (voltage), only which channel out of several.
%
% Test ID:          Integer defined in the standalone calibration data (raw format). For a specific temperature and type
% (testIdNbr)       of calibration data, it identifies a specific calibration file table with specific input/output
%                   channels (or TC-bias and bias currents), stimuli, MUX mode (+latching relay setting). Test IDs are
%                   identical for the same type of calibration data at different temperatures.
%
% RPS               Radians per second
%
%
% Types of calibration data:
% DCV (DC_VOLTAGE):         Calibration data for determining offset and slope for DC voltages.
% TF  (TRANSFER_FUNCTION):  Transfer functions (amplification and phase shift over frequencies).
%
% NOTE: The naming convention here is to use the flow of information in the physical instrument to determine what is
% "input" and "output". Therefor, for voltages, BICAS calculates the "input" signals from the "output" signals.
%
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden.
% First created 2017-12-11
%
classdef reader_DCV_TF < handle     % Explicitly declare it as a handle class to avoid IDE warnings.
% BOGIQ
% =====
% PROPOSAL: Force user to only use DCV or TF tests.
%   PROPOSAL: Internal flag that selects which.
%   PROPOSAL: Subclasses for each case.
% PROPOSAL: Implement pre-loading of data.
%   CON: Code does not know whether to load DCV or TF.
%
% PROPOSAL: Exchange calibType 'TF' --> 'FTF'.
%
% TODO: Reader for calibrating bias currents (other class).
% QUESTION: Can somehow share code with future class for bias current calibration data? Should be similar.



    properties(Access=private)
        %doPreload    = [];
        metadataList = [];
    end



    %###################################################################################################################



    methods(Access=public)

        % CONSTRUCTOR
        function obj = reader_DCV_TF(varargin)
%             if nargin > 1
%                 error('BICAS:reader_DCV_TF:IllegalArgument', 'Illegal number of arguments')
%             end
%             
%             obj.doPreload = ismember('preload', varargin);
        end


        
        function metadataList = add_test_directory(obj, cTableFilesPattern, testLogbookFile, mebTemperatureCelsius)
        % Add directory with tests (calibration tables).
        % Example: TEMP25C/4_4_DC_VOLTAGE_TEST/ or 4-5_TRANSFER_FUNCTION/.
        %
        % ARGUMENTS
        % =========
        % cTableFilesPattern    : Path as a sprintf pattern describing all test files with %02i or %03i representing
        %                         the test ID nbr. NOTE: Short test ID numbers have to be preceeded by zeroes
        %                         (due to filenaming convention).
        %                         Example: '4-5_TRANSFER_FUNCTION/SO_BIAS_AC_VOLTAGE_ID%02i_Ver_00_FS0_PAFM.txt'
        % testLogbookFile       : Path to testlogbook* file that describes and enumerates the files referred to by
        %                         cTableFilesPattern.
        % mebTemperatureCelsius : The MEB temperature at which the tests are made.
            
            testLogbookRowList = bicas.utils.read_text_file(testLogbookFile);
            metadataList = bicas.BSACT_utils.parse_testlogbook_DCV_TF(testLogbookRowList);
            
            metadataList = bicas.utils.merge_structs(metadataList, struct('mebTempCelsius', mebTemperatureCelsius));
            
            
            
            % IMPLEMENTATION NOTE: Does not read entire files since does not know whether they are DCV or TF files.
            [metadataList.filePath] = deal([]);
            for i = 1:numel(metadataList)
                filePath = sprintf(cTableFilesPattern, metadataList(i).testIdNbr);
                
                % ASSERTION: Check that files exist.
                if ~exist(filePath, 'file')
                    error('BICAS:read_text_file:Assertion', 'Can not find file "%s".', filePath)
                end
                
                metadataList(i).filePath = filePath;    % NOTE: Adding field to metadataList.
            end
            
            
            
            obj.metadataList = [obj.metadataList, metadataList];
        end
        
        
        
        function metadataList = get_metadataList(obj)
            metadataList = obj.metadataList;
        end

    end    % methods(Access=public)

    
    
    %###################################################################################################################



    methods(Static, Access=public)

        function Data = read_DCV_calib_file(filePath)
        % Read BSACT DCV file.
            
            if ~ischar(filePath); error('BICAS:reader_DCV_TF:IllegalArgument', 'Argument is not a string.'); end
            
            % Variable naming convention to distinguish columns with very similar meanings (and values):
            % EOO = Excluding output offset (offset removed)
            % IOO = Including output offset (offset not removed)
            % BST = Before stimuli
            % AST = After stimuli (after voltage drop due to stimuli)
            Data = bicas.BSACT_utils.read_BSACT_file(filePath, ...
                    {'inputBstVolt', 'outputIooVolt', 'antennaCurrentAmpere', 'inputAstVolt', 'outputEooVolt', ...
                    'hkInput1AstVolt', 'hkInput2AstVolt', 'hkInput3AstVolt'});
        end



        function Data = read_TF_calib_file(filePath)
        % Read BSACT TF file.
        
            if ~ischar(filePath) ; error('BICAS:reader_DCV_TF:IllegalArgument', 'Argument is not a string.') ; end
            
            Data = bicas.BSACT_utils.read_BSACT_file(filePath, ...
                {'freqHz', 'gainEnergyDb', 'phaseShiftDeg'});
            
            % Add fields to Data, fields which are likely to be used for plotting.
            [Data.freqRps, Data.z] = bicas.utils.convert_TF_human2math(...
                Data.freqHz, ...
                Data.gainEnergyDb, ...
                Data.phaseShiftDeg);
        end

    end    % methods(Static, Access=public)

    

    %###################################################################################################################
end
