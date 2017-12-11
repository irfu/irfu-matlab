function CTableMetadataList = parse_testlogbook(rowStrList)
%
% In raw BIAS standalone calibration data, parse subset of SO testlogbook text.
% Returns data which associates Test ID (and thus calibration data file) with various settings, e.g. mux mode.
%
% The format is approximately repeated sequences of data as below.
% """"""""
% Antenna 3, LFR Output 
% Ant 1 = GND, Ant 2 = GND, Ant 3 = Signal, Stimuli = 1Mohm
% ID68 = Mode 0 (std operation), LFR_3 = V23_DC
% ID69 = Mode 0 (std operation), LFR_5 = V23_AC, Gain = 5
% ID70 = Mode 0 (std operation), LFR_5 = V23_AC, Gain = 100
% ID71 = Mode 1 (probe 1 fails), LFR_2 = V3_DC
% """""""".
%
%
% ARGUMENTS AND RETURN VALUES
% ===========================
% CTableMetadataList : Array of structs. Each struct has fields
%   .antennaSignals   : length 3-vector. One element per antenna. Values: 0=GND (ground), 1=Signal.
%   .stimuliOhm
%   .muxMode
%   .outputChNbr
%   .inputChNbr
%   .acGain           : NaN if not explicitly stated in rowStrList.
%
% NOTE: Return result excludes latching relay (could in principle be derived from inputChNbrs sometimes).
% IMPLEMENTATION NOTE: .antennaSignals as a vector of flags i usefuls since
% (1) one can can easily check on which signals are on which inputChNbr, .e.g. all(d.antennaSignals(d.inputChNbr))
% (2) it is possible to extend the meaning of values to more alternatives than two (e.g. "Signal", "GND", "2.5 V").
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2017-10-12


CTableMetadataList = [];

i = 1;   % Row number
while i+1 <= numel(rowStrList)    % Check if reached end of file.
    % Ignore row.
    
    % Advance to next row: Header row.
    i = i+1;    
    HeaderRowSettings = parse_header_row(rowStrList{i});
    
    % Parse test rows.
    while i+1 <= numel(rowStrList)    % Check if reached end of file.
        % Advance to next row.
        i = i+1;
        
        rowStr = rowStrList{i};
        RowCTableMetadata = parse_test_row(rowStr);
        if isempty(RowCTableMetadata)
            break
        end
        CTableMetadataList = [CTableMetadataList, ...
            bicas.utils.merge_structs(HeaderRowSettings, RowCTableMetadata)];
        %CTableMetadataList(end+1) = bicas.utils.merge_structs(HeaderRowSettings, RowCTableMetadata);
    end
    
end

end



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
% RETURN VALUE
% ============
% CTableMetadata : Struct with derived values. Empty if not the intended type of row.
%
% Example rowStr that function should be able to handle.
% "ID00 = Mode 0 (std operation), LFR_1 = V1_DC"
% "ID01 = Mode 0 (std operation), LFR_2 = V12_DC*"
% "ID02 = Mode 0 (std operation), LFR_4 = V12_AC*, Gain = 5"
% "ID43 = Mode 4 (cal mode 0), TDS_1 = V1_DC"

% DEBUG
%fprintf('rowStr = "%s"\n', rowStr);

% IMPLEMENTATION NOTE: Must use %d, not %i which can be interpreted as octal. Therefore, use 'ID%d'.
% """"%i    Base determined from the values. Defaults to base 10. If initial digits are 0x or 0X, it is base 16. If initial digit
% is 0, it is base 8.""""
% Ex: sscanf('ID07', 'ID%i') ==> 7
%     sscanf('ID08', 'ID%i') ==> 0
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
