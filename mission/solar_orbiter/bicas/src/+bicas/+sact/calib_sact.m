classdef calib_sact < handle     % Explicitly declare it as a handle class to avoid IDE warnings.
%
% IMPLEMENTATION NOTE:
% UNFINISHED CODE. NOT USED BY MAIN PROGRAM YET. NOT ENTIRELY CLEAR WHAT IT SHOULD DO.
%
% Class for code and data related to classifying, locating, and reading raw BIAS standalone calibration data files
% in the format, filenames, and directory structures they have from the calibration campaign.
%
% The non-trivial parts are centered on deriving metadata for every calibration table
%
% IMPLEMENTATION NOTE: The way of specifying files and parameters describing (classifyin) the calibration files has
% deliberately been kept as close to the format of the original files as possible.
%
% IMPLEMENTATION NOTE: Loading all calibration files rather than just their fits, or the relevant subsets could be
% useful:
% (1) Can be used for implementing more specific uses of calibration data
%       Example 1: Select AC (spectrum) transfer functions (TF) based on combination of "channel" and high gain/low
%       gain, rather than just one for AC low gain and AC high gain as expected.
%       Example 2: Interpolate table instead of linear function (offset & slope; y=k*x+m).
%       Example 3: Derive fits (DC) from calibration tables directly rather than using separately tabulated offsets and
%       slopes.
% (2) Manually analyzing the calibration data itself (plotting; making one's own fits using other table columns; comparing
% channels, checking temperature dependence).
% (3) Converting the calibration files to more suitable formats (e.g. XML).
%
%
% VARIABLE NAMING CONVENTIONS / DEFINITIONS OF TERMS
% ==================================================
% SACT              StandAlone Calibration Tables
%
% CTable:           Calibration table. Equal to the contents of one BIAS standalone calibration file.
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
% Calibration type: The type of calibration data. One of the following strings
% (calibType)       'DCV' (DC_VOLTAGE):         Calibration data for determining offset and slope for DC voltages.
%                   'TF'  (TRANSFER_FUNCTION):  Transfer functions (amplification and phase shift over frequencies).
%
% NOTE: Not always obvious what counts as "input" and "output". The convention here is to use the flow of information
% in the physical instrument. Therefor, for voltages, BICAS calculates the "input" signals from the "output" signals.
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2017-03-15


% BOGIQ
% =====
% PROPOSAL: Change to better name than SACT.
%   PROPOSAL: BSACT  = BIAS StandAlone Calibration Tables (BSCT?)
%   PROPOSAL: RBSACT = Raw BIAS StandAlone Calibration Tables (RBSCT?)
%   PROPOSAL: RBSCD  = Raw BIAS Standalone Calibration Data
%   PROPOSAL: RBSACTR, RBSCTR, RBSCR = Raw BIAS Stand-Alone Calibration Tables Reader
%   PROPOSAL: Terms: Raw, BIAS, access/read/reader, standalone calibration
% PROPOSAL: Change name of class/file. (XXXX=acronym)
%   PROPOSAL: XXXX_class, XXXX_object, XXXX_reader, XXXXX (e.g. rbscr).
%
% PROPOSAL: Change name inputChNbrs --> probeNbrs
%   CON: We want analogous signal names ==> "probeSignalVolt" : Not good name for both single and diff?
%       PROPOSAL: probeSingleSignalVolt
%                 probeDiffSignalVolt
%
% PROPOSAL: Move SACT out of BICAS?
%   PRO: BICAS will not use SACT itself for calibration in the processing.
%       CON: It might be used for deriving calibration files (via fits of standalone calibration data).
%   CON: It is a tool related to BICAS.
%       PRO: Should be properly stored and archived somewhere.
%   CON: SACT uses BICAS code: bicas.calibration.*, bicas.utils.*.
%
% PROPOSAL: Calibration type/calibType effectively abolished. Remove.
%   CON: Redo into non-variable. Naming convention, context needed.
%
% PROPOSAL: Specify file names with globbing, regexp.
% PROPOSAL: Move initialization functions into constructor.
% PROPOSAL: Merge DcvCTableMetadataList, TfCTableMetadataList (add extra field to distinguish). Similar to
%           FilePathPatternList.
%   PRO: There is identical code that works on both.
%
% PROPOSAL: Shorten Degrees-->Deg, RadPerSec-->Rps, frequency-->freq
%
% PROPOSAL: Exchange calibType 'TF' --> 'FTF'.

    properties(Access=private)

        standaloneCalibFilesRootPath = [];

        % Struct array.
        % List of file path patterns for ALL calibration data types.
        FilePathPatternList = [];
        
        DcvCTableMetadataList  = [];
        TfCTableMetadataList   = [];
        
    end
    
    %###################################################################################################################

    methods(Access=public)
    
        % CONSTRUCTOR
        function obj = calib_sact(standaloneCalibFilesRootPath, dcvRawTestLogbookPath, tfRawTestLogbookPath)
            % ARGUMENTS
            % =========
            % standaloneCalibFilesRootPath : Path to the root directory of BIAS standalone calibration files
            % dcvRawTestLogbookPath        : Path to extract of the test logbook for "4_4_DC_VOLTAGE_TEST/" tests in BIAS
            %                                standalone calibration files.
            % tfRawTestLogbookPath         : Path to extract of the test logbook for "4_5_TRANSFER_FUNCTION/"
            %                                (AC frequency transfer functions) tests in BIAS standalone calibration files.
            % NOTE: "extract of the test logbook" here refers to the section of the test logbook files that describes
            % which test corresponds to a given test ID. Only these sections can presently be automatically parsed.
            
            obj.standaloneCalibFilesRootPath = standaloneCalibFilesRootPath;            
            
            obj.FilePathPatternList = bicas.sact.calib_sact.init_FilePathPatternLists();
            
            [obj.DcvCTableMetadataList, obj.TfCTableMetadataList] = bicas.sact.calib_sact.init_cTableMetadataList(...
                dcvRawTestLogbookPath, tfRawTestLogbookPath);
        end



        function filePath = get_SACT_table_file_path(obj, calibType, mebTempCelsius, testIdNbr)

            % Select file pattern record to use.
            tempList = bicas.utils.select_array_structs(obj.FilePathPatternList, 'calibType',      {calibType});
            tempList = bicas.utils.select_array_structs(tempList,                'mebTempCelsius', {mebTempCelsius});
            
            % ASSERTION
            if numel(tempList) ~= 1
                error('BICAS:get_SACT_table_file_path:Assertion', 'Does not find exactly one path. numel(tempList)=%i', numel(tempList))
            end
            
            relativeFilePath = sprintf(tempList.filePathPattern, testIdNbr);
            filePath = fullfile(obj.standaloneCalibFilesRootPath, relativeFilePath);
        end
        
        
        
        function List = get_DcvCTableMetadataList(obj)
            List = obj.DcvCTableMetadataList;
        end

        function List = get_TfCTableMetadataList(obj)
            List = obj.TfCTableMetadataList;
        end

        
        
%         function calibTable = get_DC_VOLTAGE_calib_table(obj, varargin)
%         % Obtain a specific DC_VOLTAGE calibration table.
%         %
%         % ASSUMES: Exact matches for stimuli and temperature. No interpolation.
%         
%             if length(varargin) == 4        
%                 [calibTable, iTable] = bicas.utils.select_array_structs(...
%                     obj.calibTablesList, ...
%                     'mebTempCelsius', {varargin{1}}, ...
%                     'stimuliOhm',     {varargin{2}}, ...
%                     'inputChNbr',     {varargin{3}}, ...
%                     'outputChNbr',    {varargin{4}});
%             elseif length(varargin) == 2
%                 % Case is useful for plotting calibration tables, and if one wants to iterate over testIdNbr.
%                 [calibTable, iTable] = bicas.utils.select_array_structs(...
%                     obj.calibTablesList, ...
%                     'mebTempCelsius', {varargin{1}}, ...
%                     'testIdNbr',      {varargin{2}});
%             else
%                 error('BICAS:get_DC_VOLTAGE_calib_table:Assertion', 'Wrong number of arguments.')
%             end
%             
%             if numel(calibTable) ~= 1
%                 error('BICAS:get_DC_VOLTAGE_calib_table:Assertion', 'Can not find exactly one calibration table for the given values.')
%             end
%             
%             if isempty(calibTable.data)
%                 filePath = obj.get_SACT_table_file_path('DCV', calibTable.mebTempCelsius, calibTable.testIdNbr);
%                 
%                 obj.calibTablesList(iTable).data = read_DC_VOLTAGE_calib_file(filePath);
%                 
%                 calibTable = obj.calibTablesList(iTable);
%             end            
%         end

    end    % methods(Access=public)

    %###################################################################################################################

    methods(Static, Access=public)

        function Data = read_DC_VOLTAGE_calib_file(filePath)
        % Read BIAS standalone calibration file for DC voltage offset+slope.
        %
        % "DC_VOLTAGE" refers to the 4_4_DC_VOLTAGE_TEST// directories in the raw BIAS standalone calibration
        % data.
            
            % Variable naming convention to distinguish columns with very similar meanings (and values):
            % EOO = Excluding output offset (offset removed)
            % IOO = Including output offset (offset not removed)
            % BST = Before stimuli
            % AST = After stimuli (after voltage drop due to stimuli)
            Data = bicas.sact.calib_sact.read_sact_file(filePath, ...
                    {'inputBstVolt', 'outputIooVolt', 'antennaCurrentAmpere', 'inputAstVolt', 'outputEooVolt', ...
                    'hkInput1AstVolt', 'hkInput2AstVolt', 'hkInput3AstVolt'});
        end



        function Data = read_TRANSFER_FUNCTION_file(filePath)
        % Read BIAS standalone calibration file describing a transfer function (over frequency).
        % 
        % "TRANSFER_FUNCTION" refers to the 4_5_TRANSFER_FUNCTION/ directories in the raw BIAS standalone calibration
        % data.
        
        % TODO: Assertions? Implicit check on file format.
        
            Data = bicas.sact.calib_sact.read_sact_file(filePath, ...
                {'frequencyHz', 'gainEnergyDb', 'phaseShiftDegrees'});
            
            %UPPER_WRAP_BOUNDARY_DEG = 10;
            %Data.phaseShiftDegrees = wrapTo360(Data.phaseShiftDegrees + 360 - UPPER_WRAP_BOUNDARY_DEG) + UPPER_WRAP_BOUNDARY_DEG;            
            %wrapTo180(Data.phaseShiftDegrees(0) / 2)
            %Data.phaseShiftDegrees = Data.phaseShiftDegrees
            
            % Add fields to Data.
            [Data.frequencyRadPerSec, Data.z] = bicas.calibration.convert_TF_Hz_energy_dB_phase(...
                Data.frequencyHz, ...
                Data.gainEnergyDb, ...
                Data.phaseShiftDegrees);
        end

    end   % methods(Static, Access=public)

    %###################################################################################################################

    methods(Static, Access=private)
        
        
        
        function FilePathPatternList = init_FilePathPatternLists
        % Initialize path patterns that can be used to derive the paths of specific standalone calibration files (only
        % the test ID is missing).
        
            % PROPOSAL: Load from BICAS settings or constructor?
            temp = [];
            temp = [temp, set_record('DCV', -25, fullfile('TEMP-25C', '4_4_DC_VOLTAGE_TEST', 'SO_BIAS_DC_VOLTAGE_ID%i_Ver_00_FM1_-25_4.4.txt'))];
            temp = [temp, set_record('DCV',   0, fullfile('TEMP0C',   '4_4_DC_VOLTAGE_TEST', 'SO_BIAS_DC_VOLTAGE_ID%i_Ver_00_FM1_0_4.4.txt'))];
            temp = [temp, set_record('DCV',  25, fullfile('TEMP25C',  '4_4_DC_VOLTAGE_TEST', 'SO_BIAS_DC_VOLTAGE_ID%i_Ver_00_FM1_+25C.txt'))];
            temp = [temp, set_record('DCV',  50, fullfile('TEMP50C',  '4_4_DC_VOLTAGE_TEST', 'SO_BIAS_DC_VOLTAGE_ID%i_Ver_00_FM1_+50c.txt'))];
            temp = [temp, set_record('DCV',  70, fullfile('TEMP70C',  '4_4_DC_VOLTAGE_TEST', 'SO_BIAS_DC_VOLTAGE_ID%i_Ver_00_FM1.txt'))];

            temp = [temp, set_record('TF', -25, fullfile('TEMP-25C', '4_5_TRANSFER_FUNCTION', 'SO_BIAS_AC_VOLTAGE_ID%02i_Ver_00_FM1_-25_4.5.txt'))];
            temp = [temp, set_record('TF',   0, fullfile('TEMP0C',   '4_5_TRANSFER_FUNCTION', 'SO_BIAS_AC_VOLTAGE_ID%02i_Ver_00_FM1_0_4.5.txt'))];
            temp = [temp, set_record('TF',  25, fullfile('TEMP25C',  '4_5_TRANSFER_FUNCTION', 'SO_BIAS_AC_VOLTAGE_ID%02i_Ver_00_FM1_+25C.txt'))];
            temp = [temp, set_record('TF',  50, fullfile('TEMP50C',  '4_5_TRANSFER_FUNCTION', 'SO_BIAS_AC_VOLTAGE_ID%02i_Ver_00_FM1_+50c.txt'))];
            temp = [temp, set_record('TF',  70, fullfile('TEMP70C',  '4_5_TRANSFER_FUNCTION', 'SO_BIAS_AC_VOLTAGE_ID%02i_Ver_00_FM1.txt'))];
            FilePathPatternList = temp;

            % NESTED FUNCTION
            function Record = set_record(calibType, mebTempCelsius, filePathPattern)
                Record = struct(...
                    'calibType',       calibType, ...
                    'mebTempCelsius',  mebTempCelsius, ...
                    'filePathPattern', filePathPattern);
            end
        end
        
        

        function [DcvCTableMetadataList, TfCTableMetadataList] = init_cTableMetadataList(...
                dcvRawTestLogbookPath, tfRawTestLogbookPath)
        % Function for generating list of records representing all calibration file tables and their metadata.
        % Presently only covers 4_4_DC_VOLTAGE_TEST, 4_5_TRANSFER_FUNCTION.

            mebTempCelsiusList = [-25, 0, 25, 50, 70];   % Should ideally be derived from obj.FilePathPatternList.
            
            dcvRawTestLogbookRowList = bicas.utils.read_text_file(dcvRawTestLogbookPath);
            tfRawTestLogbookRowList  = bicas.utils.read_text_file(tfRawTestLogbookPath);
            
            dcvTestLogbookList = bicas.sact.parse_testlogbook(dcvRawTestLogbookRowList);
            tfTestLogbookList  = bicas.sact.parse_testlogbook(tfRawTestLogbookRowList);
            
            DcvCTableMetadataList = [];
            TfCTableMetadataList  = [];
            for mebTempCelsius = mebTempCelsiusList
                DcvCTableMetadataList = [DcvCTableMetadataList, ...
                    bicas.utils.merge_structs(dcvTestLogbookList, struct('mebTempCelsius', mebTempCelsius))];
                
                TfCTableMetadataList = [TfCTableMetadataList, ...
                    bicas.utils.merge_structs(tfTestLogbookList, struct('mebTempCelsius', mebTempCelsius))];
            end
            
            DcvCTableMetadataList = bicas.sact.calib_sact.derive_extra_cTable_metadata(DcvCTableMetadataList);
            TfCTableMetadataList  = bicas.sact.calib_sact.derive_extra_cTable_metadata(TfCTableMetadataList);

        end


        
        function data = read_sact_file(filePath, columnFieldNamesList)
        % Read raw BIAS standalone calibration file.
        %
        % Reads a file file on the format of the files generated by the actual BIAS standalone calibration campaign, and
        % NOT files later converted into any more "convenient" format. It effectively reads a numeric table in
        % a text file with a certain comment style and labels the columns according to arguments.
        %
        %
        % ARGUMENTS AND RETURN VALUE
        % ==========================
        % columnFieldNamesList : Cell array of strings, one per file table column (left to right).
        % data                 : struct with one field per column.
        %
        %
        % IMPLEMENTATION NOTE: Code needs to be able to handle files with different numbers of columns.
        
            fId = fopen(filePath, 'r');
            if fId == -1
                error('BICAS:read_sact_file:PathNotFound', 'Can not open file "%s".', filePath)
            end
            fileContents = textscan(fId, '%f', 'CommentStyle', 'mheader');   % Format specification will read all numbers into one long 1D array.
            fileContents = fileContents{1};
            fclose(fId);

            nColumns = length(columnFieldNamesList);
            
            % ASSERTION. Tries to check for the number of columns, but is not a perfect test.
            if mod(length(fileContents), nColumns) ~= 0
                error('BICAS:read_sact_file:UnexpectedFileFormat', 'Number of specificed columns does not match file contents.')
            end
            fileContents = reshape(fileContents, [nColumns, length(fileContents)/nColumns])';
            
            data = struct;
            for iColumn = 1:nColumns
                fieldName = columnFieldNamesList{iColumn};
                data.(fieldName) = fileContents(:, iColumn);
            end
        end



        function CTableMetadata = derive_extra_cTable_metadata(CTableMetadata)
        % ARGUMENTS
        % =========
        % CTableMetadata : Calibration table metadata struct array.
        
        % PROPOSAL: Add latchingRelay, isDiff.
        % PROPOSAL: Add phaseShiftNormalizedHz.
        %   NOTE: Only exists for TF.
        
            % Create empty new fields.
            [CTableMetadata.invertedInput]   = deal([]);
            [CTableMetadata.commonModeInput] = deal([]);
        
            for i = 1:numel(CTableMetadata)
                isDiff         = numel(CTableMetadata(i).inputChNbr) == 2;
                inputChSignals = CTableMetadata(i).antennaSignals(CTableMetadata(i).inputChNbr);
                
                % True iff input is inverted due to choice of which antenna a diff is taken.
                CTableMetadata(i).invertedInput   = isDiff && all(inputChSignals == [0,1]);
                % True iff diff and signal on both antennas.
                CTableMetadata(i).commonModeInput = isDiff && all(inputChSignals == [1,1]);
            end            
        end
        
        
        
        % PROPOSAL: Move out of calib_sact. Separate function file?
%         function MuxSettings = get_mux_settings(muxMode, latchingRelay)
%         % EXPERIMENTAL
%         % Function which describes the demultiplexer conceptually as in Table 4, RPW-SYS-MEB-BIA_SPC-00001-IRF,
%         % "Specifications of RPW/BIAS".
%         %
%         % ARGUMENTS AND RETURN VALUES
%         % ===========================
%         % muxMode
%         % latchingRelay : Latching relay flag as mentioned in Table 4, RPW-SYS-MEB-BIA_SPC-00001-IRF. String.
%         %                   '12' = Use V12_DC/AC.
%         %                   '13' = Use V13_DC/AC.
%         % MuxSettings   : struct with which input channels are used for which output channels.
%         %                 Special case: empty : There is no regular signal (2.5 V, GND).
%         %   .BIAS_1_DC_inputChNbrs :
%         %   .BIAS_2_DC_inputChNbrs :
%         %   .BIAS_3_DC_inputChNbrs :
%         %   .BIAS_4_AC_inputChNbrs :
%         %   .BIAS_5_AC_inputChNbrs :
%         %
%         % It is conceiveable that this function should be moved out of here and that
%         % dm_processing_functions.simple_demultiplex_subsequence should be based on it.
%         %
%         % TODO: Latching relay
%         
%             switch(latchingRelay)
%                 case '12'
%                     BIAS_DCAC_inputChNbrs_latchingRelay = [1,2];
%                 case '13'
%                     BIAS_DCAC_inputChNbrs_latchingRelay = [1,3];
%                 otherwise
%             end
% 
%             switch(muxMode)
%                 case 0   % "Standard operation" : We have all information.
%                     
%                     % Summarize the INPUT DATA we have.
%                     MuxSettings.BIAS_1_DC_inputChNbrs = [1];
%                     MuxSettings.BIAS_2_DC_inputChNbrs = BIAS_DCAC_inputChNbrs_latchingRelay;
%                     MuxSettings.BIAS_3_DC_inputChNbrs = [2 3];
%                     MuxSettings.BIAS_4_AC_inputChNbrs = BIAS_DCAC_inputChNbrs_latchingRelay;
%                     MuxSettings.BIAS_5_AC_inputChNbrs = [2 3];
%                     
%                 case 1   % Probe 1 fails
%                     
%                     MuxSettings.BIAS_1_DC_inputChNbrs = [2];
%                     MuxSettings.BIAS_2_DC_inputChNbrs = [3];
%                     MuxSettings.BIAS_3_DC_inputChNbrs = [2 3];
%                     MuxSettings.BIAS_4_AC_inputChNbrs = BIAS_DCAC_inputChNbrs_latchingRelay;
%                     MuxSettings.BIAS_5_AC_inputChNbrs = [2 3];
%                     % Input.BIAS_4 unavailable.
%                     
%                 case 2   % Probe 2 fails
%                     
%                     MuxSettings.BIAS_1_DC_inputChNbrs = [1];
%                     MuxSettings.BIAS_2_DC_inputChNbrs = [3];
%                     MuxSettings.BIAS_3_DC_inputChNbrs = BIAS_DCAC_inputChNbrs_latchingRelay;
%                     MuxSettings.BIAS_4_AC_inputChNbrs = BIAS_DCAC_inputChNbrs_latchingRelay;
%                     MuxSettings.BIAS_5_AC_inputChNbrs = [2 3];
%                     % Input.BIAS_5 unavailable.
%                     
%                 case 3   % Probe 3 fails
%                     
%                     MuxSettings.BIAS_1_DC_inputChNbrs = [1];
%                     MuxSettings.BIAS_2_DC_inputChNbrs = [2];
%                     MuxSettings.BIAS_3_DC_inputChNbrs = BIAS_DCAC_inputChNbrs_latchingRelay;
%                     MuxSettings.BIAS_4_AC_inputChNbrs = BIAS_DCAC_inputChNbrs_latchingRelay;
%                     MuxSettings.BIAS_5_AC_inputChNbrs = [2 3];
%                     % Input.BIAS_5 unavailable.
%                     
%                 case 4   % Calibration mode 0
%                     
%                     MuxSettings.BIAS_1_DC_inputChNbrs = [1];
%                     MuxSettings.BIAS_2_DC_inputChNbrs = [2];
%                     MuxSettings.BIAS_3_DC_inputChNbrs = [3];
%                     MuxSettings.BIAS_4_AC_inputChNbrs = BIAS_DCAC_inputChNbrs_latchingRelay;
%                     MuxSettings.BIAS_5_AC_inputChNbrs = [2 3];
% 
%                 case {5,6,7}   % Calibration mode 1/2/3
%                     
%                     MuxSettings.BIAS_1_DC_inputChNbrs = [];
%                     MuxSettings.BIAS_2_DC_inputChNbrs = [];
%                     MuxSettings.BIAS_3_DC_inputChNbrs = [];
%                     MuxSettings.BIAS_4_AC_inputChNbrs = BIAS_DCAC_inputChNbrs_latchingRelay;
%                     MuxSettings.BIAS_5_AC_inputChNbrs = [2 3];
%                     
%                 otherwise
%                     if isnan(muxMode)
%                         ;   % Do nothing. Allow the default values (NaN) to be returned.
%                     else
%                         error('BICAS:data_manager:Assertion:IllegalArgument:DatasetFormat', 'Illegal argument value for mux_set.')
%                     end
%             end    % switch
%         end    % function
        
    end    % methods(Static, Access=private)
        
end

