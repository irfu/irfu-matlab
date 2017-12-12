function metadataList = parse_testlogbook_DCV_TF(rowStrList)
%
% In the raw BIAS standalone calibration data, parse testlogbook file for
% type DC_VOLTAGE and TRANSFER_FUNCTION to obtain metadata for calibration table files.
% Returns data which associates Test ID (and thus calibration data file) with various settings, e.g. mux mode.
%
%
% TEXT FILE FORMAT
% ================
% The format is a text file format in the following sequence.
% (1) Various human-readable text (ignored),
% (2) Repeated sequences of data as below
% """"""""
% Antenna 3, LFR Output 
% Ant 1 = GND, Ant 2 = GND, Ant 3 = Signal, Stimuli = 1Mohm
% ID68 = Mode 0 (std operation), LFR_3 = V23_DC
% ID69 = Mode 0 (std operation), LFR_5 = V23_AC, Gain = 5
% ID70 = Mode 0 (std operation), LFR_5 = V23_AC, Gain = 100
% ID71 = Mode 1 (probe 1 fails), LFR_2 = V3_DC
% """""""".
% (3) Various human-readable text and plotting and fitting log messages (ignored).
%
%
%
% ARGUMENTS AND RETURN VALUES
% ===========================
% rowList             : Cell array of strings representing rows in text file.
% metadataList        : Array of structs. Each struct has fields:
%   .testIdNbr
%   .antennaSignals   : Length 3-vector. [i] = Value for antenna i. Values: 0=GND (ground), 1=Signal.
%   .stimuliOhm
%   .muxMode
%   .outputChNbr
%   .inputChNbr       : Scalar array length 1 or 2. Contains antenna number(s). inputChNbr(1) < inputChNbr(2).
%   .acGain           : NaN if not explicitly stated in rowStrList.
%   .invertedInput    : True iff diff AND inputChNbr(1)=GND and inputChNbr(2)=Signal.
%   .commonModeInput  : True iff diff AND common mode (2 signals).
%   (but not .mebTempCelsius, .filePath)
%
%
%
% NOTES
% =====
% NOTE: Return result excludes latching relay (could in principle be derived from inputChNbrs sometimes).
% NOTE: The class name is chosen to reflect the types of calibration data that it may contain in anticipation of 
% eventually creating another analogous class for other calibration data (bias current calibration data).
% IMPLEMENTATION NOTE: .antennaSignals as a vector of flags is useful since
% (1) one can can easily check on which signals are on which inputChNbr, .e.g. all(d.antennaSignals(d.inputChNbr))
% (2) it is possible to extend the meaning of values to more alternatives than two (e.g. "Signal", "GND", "2.5 V").
% IMPLEMENTATION NOTE: Takes list of rows as argument instead of file path, partly to make automated testing easier.
%
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2017-10-12

%==============================================================================================
% BOGIQ
% -----
% PROPOSAL: Rename return variable to imply that fields will later be added to it.
% PROPOSAL: Check for existence of Header1 rows, just for confirmation.
% PROPOSAL: inputChNbr --> inputChNbrs (plural)
% PROPOSAL: commonMode --> isCommonMode
% PROPOSAL: invertedInput --> hasInvertedInput
% PROPOSAL: Should become part of reader_DCV_TF?
%==============================================================================================


% Strings representing different types of rows in testlogbook file:
% =================================================================
% NoDataPrefix : Any row before the first "Header2" row.
% NoDataSuffix : Any row after  the last  "Test"    row.
% Header1      : E.g. "Antenna 3, LFR Output"           % (skipped by algorithm by incrementing row)
% Header2      : E.g. "Ant 1 = GND, Ant 2 = GND, Ant 3 = Signal, Stimuli = 1Mohm"
% Test         : E.g. "ID68 = Mode 0 (std operation), LFR_3 = V23_DC"

metadataList = [];
iRow = 1;
expectedRowType = 'NoDataPrefix_Header2';   % NOTE: Underscore should be interpreted as "or".
while iRow <= numel(rowStrList)    % Check if reached end of file.
    % CASE: Row iRow exists.
    
    rowStr = rowStrList{iRow};
    
    %fprintf('expectedRowType = %s\n', expectedRowType);   % DEBUG
    %fprintf('rowStr = "%s"\n', rowStr);
    
    switch expectedRowType

        case 'NoDataPrefix_Header2'
            expectedRowType = 'Test';
            try
                HeaderRowSettings = parse_header_row(rowStrList{iRow});
            catch
                expectedRowType = 'NoDataPrefix_Header2';
            end
            iRow = iRow + 1;
            
        case 'NoDataSuffix_Header2'
            expectedRowType = 'Test';
            try
                HeaderRowSettings = parse_header_row(rowStrList{iRow});
            catch
                expectedRowType = 'NoDataSuffix';
            end
            iRow = iRow + 1;
            
        case 'Test'
            
            RowCTableMetadata = parse_test_row(rowStr);
            if isempty(RowCTableMetadata)
                expectedRowType = 'NoDataSuffix_Header2';
                iRow = iRow + 1;   % Skip row since skipped row should have been a Header1 row.
            else
                metadataList = [...
                    metadataList, ...
                    bicas.utils.merge_structs(HeaderRowSettings, RowCTableMetadata)];   % Works with metadataList==[];
                expectedRowType = 'Test';
                iRow = iRow + 1;
            end
            
        case 'NoDataSuffix'
            break

        otherwise
            error('BICAS:parse_testlogbook', 'Assertion');

    end    % switch
    
end    % while



metadataList = derive_extra_cTable_metadata(metadataList);

end % functions


function Settings = parse_header_row(rowStr)
% Parse testlogbook "header row", e.g. "Ant 1 = GND, Ant 2 = GND, Ant 3 = Signal, Stimuli = 1Mohm".

    Settings.antennaSignals = [...
        map_regex_to_values(rowStr, 'Ant 1 = GND', 0, 'Ant 1 = Signal', 1), ...
        map_regex_to_values(rowStr, 'Ant 2 = GND', 0, 'Ant 2 = Signal', 1), ...
        map_regex_to_values(rowStr, 'Ant 3 = GND', 0, 'Ant 3 = Signal', 1)];
    
    Settings.stimuliOhm = map_regex_to_values(rowStr, '100kohm', 1e5, '1Mohm', 1e6);
end



function CTableMetadata = parse_test_row(rowStr)
% Interpret a row with information about a specific test, e.g. "ID68 = Mode 0 (std operation), LFR_3 = V23_DC".
%
%
% RETURN VALUE
% ============
% CTableMetadata : Struct with derived values. Empty if not the intended type of row.
%
%
% Example rowStr that function should be able to handle.
% "ID00 = Mode 0 (std operation), LFR_1 = V1_DC"
% "ID01 = Mode 0 (std operation), LFR_2 = V12_DC*"
% "ID02 = Mode 0 (std operation), LFR_4 = V12_AC*, Gain = 5"
% "ID43 = Mode 4 (cal mode 0), TDS_1 = V1_DC"


% DEBUG
%fprintf('rowStr = "%s"\n', rowStr);

%================================================================================================================
% IMPLEMENTATION NOTE: Must use %d, not %i which can be interpreted as octal. Therefore, use 'ID%d'.
% """"%i    Base determined from the values. Defaults to base 10. If initial digits are 0x or 0X, it is base 16.
% If initial digit is 0, it is base 8.""""
% Ex: sscanf('ID07', 'ID%i') ==> 7
%     sscanf('ID08', 'ID%i') ==> 0
%================================================================================================================
CTableMetadata.testIdNbr = find_parse_nbr(rowStr, 'ID[0-9]*',      'ID%d', 1);
if isnan(CTableMetadata.testIdNbr)
    % CASE: (Assumption) This row is not a test settings row.
    CTableMetadata = [];
    return 
end

CTableMetadata.muxMode     = find_parse_nbr(rowStr, 'Mode [0-7]*',     'Mode %d',   0);
CTableMetadata.outputChNbr = map_regex_to_values(rowStr, ...
    'LFR_1', 1, ...
    'LFR_2', 2, ...
    'LFR_3', 3, ...
    'LFR_4', 4, ...
    'LFR_5', 5, ...
    'TDS_1', 1, ...
    'TDS_2', 2, ...
    'TDS_3', 3 ...
);
CTableMetadata.acGain      = find_parse_nbr(rowStr, 'Gain = [015]*',   'Gain = %d', 1);
CTableMetadata.inputChNbr = map_regex_to_values(rowStr, ...
    'V1_', 1, ...
    'V2_', 2, ...
    'V3_', 3, ...
    'V12_', [1 2], ...
    'V13_', [1 3], ...
    'V23_', [2 3] ...
);

% DEBUG
%fprintf('CTableMetadata.testIdNbr = "%i"\n', CTableMetadata.testIdNbr);
end



function x = map_regex_to_values(str, varargin)
% Given a string and a list of pairs (regex patterns, value), return the value for the regex pattern which is
% actually contained in the string.
%
% ARGUMENTS
% =========
% varargin : pairs of arguments: (regex pattern) + (value).
%
% NOTE: No assertion on multiple matches. First match is used.

i = 1;
while i <= numel(varargin)
    regexPattern = varargin{i};
    value = varargin{i+1};
    
    if regexp(str, regexPattern, 'start')
        x = value;
        return
    end
    i = i + 2;
end

error('BICAS:parse_testlogbook:Assertion', 'Can not find regex in string.')
end



function x = find_parse_nbr(str, regexPattern, regexMatchSscanfFormat, canBeNonExistent)
% Search string for a regex. The string that matches the regex is then parsed with sscanf.
%
% RETURN VALUE
% ============
% x : The numeric value in regexMatchSscanfFormat. NaN if there was no regex match and canBeNonExistent==true.

regexStrMatch = regexp(str, regexPattern, 'match');

if isempty(regexStrMatch)
    if canBeNonExistent
        x = NaN;
    else
        error('BICAS:find_parse_nbr', 'Can not find regexPattern="%s" in str="%s".', regexPattern, str)
    end
else
    x = sscanf(regexStrMatch{1}, regexMatchSscanfFormat);
    if isempty(x)
        if canBeNonExistent
            error('BICAS:find_parse_nbr', ...
                'sscanf can not interpret regexStrMatch{1}="%s" as regexMatchSscanfFormat="%s".', ...
                regexStrMatch{1}, regexMatchSscanfFormat)
        else
            x = NaN;
        end
    end
end
end



function metadataList = derive_extra_cTable_metadata(metadataList)
% Derive and add certain extra fields and add to metadataList.
%
% The fields which are added are:
%   invertedInput   :
%   commonModeInput :
%
%
% ARGUMENTS
% =========
% metadataList : Calibration table metadata struct array.

% PROPOSAL: Add latchingRelay, isDiff.


% Create empty new fields.
[metadataList.invertedInput]   = deal([]);
[metadataList.commonModeInput] = deal([]);

for i = 1:numel(metadataList)
    isDiff         = numel(metadataList(i).inputChNbr) == 2;
    inputChSignals = metadataList(i).antennaSignals(metadataList(i).inputChNbr);

    % True iff diff and input is inverted due to choice of which antenna a diff is taken.
    metadataList(i).invertedInput   = isDiff && all(inputChSignals == [0,1]);   % [GND, Signal]

    % True iff diff and signal on both antennas.
    % NOTE: Relies on inputChNbr(1) < inputChNbr(2).
    metadataList(i).commonModeInput = isDiff && all(inputChSignals == [1,1]);   % [Signal, Signal]
end
end
