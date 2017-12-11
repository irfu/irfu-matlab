function parse_testlogbook_TEST
%
% Automatic test code for function parse_testlogbook
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2017-10-12


exp = struct(...
    'antennaSignals',  {}, ...
    'stimuliOhm',      {}, ...
    'testIdNbr',       {}, ...
    'muxMode',         {}, ...
    'outputChNbr',     {}, ...
    'inputChNbr',      {}, ...
    'acGain',          {});

rowList = {};

rowList = [rowList, {...
'Antenna 1-2 Common Mode, LFR Output', ...
'Ant 1 = Signal, Ant 2 = Signal, Ant 3 = GND, Stimuli = 100kohm', ...
'ID09 = Mode 0 (std operation), LFR_2 = V12_DC*', ...
'ID33 = Mode 0 (std operation), LFR_4 = V12_AC*, Gain = 5', ...
'ID34 = Mode 0 (std operation), LFR_4 = V12_AC*, Gain = 100', ...
'ID35 = Mode 3 (probe 3 fails), LFR_3 = V12_DC*' }];
exp(end+1) = create_calib_table_record( 9, 0, 2, [1 2]);
exp(end+1) = create_calib_table_record(33, 0, 4, [1 2], 5);
exp(end+1) = create_calib_table_record(34, 0, 4, [1 2], 100);
exp(end+1) = create_calib_table_record(35, 3, 3, [1 2]);
exp = set_field(exp, 4, 'antennaSignals', [1 1 0]);
exp = set_field(exp, 4, 'stimuliOhm', 1e5);

rowList = [rowList, {...
'Antenna 2-3 Common Mode, LFR Output ', ...
'Ant 1 = GND, Ant 2 = Signal, Ant 3 = Signal, Stimuli = 1Mohm', ...
'ID36 = Mode 0 (std operation), LFR_3 = V23_DC', ...
'ID37 = Mode 0 (std operation), LFR_5 = V23_AC, Gain = 5', ...
'ID38 = Mode 0 (std operation), LFR_5 = V23_AC, Gain = 100', ...
'ID39 = Mode 1 (probe 1 fails), LFR_3 = V23_DC' }];
exp(end+1) = create_calib_table_record(36, 0, 3, [2 3]);
exp(end+1) = create_calib_table_record(37, 0, 5, [2 3], 5);
exp(end+1) = create_calib_table_record(38, 0, 5, [2 3], 100);
exp(end+1) = create_calib_table_record(39, 1, 3, [2 3]);
exp = set_field(exp, 4, 'antennaSignals', [0 1 1]);
exp = set_field(exp, 4, 'stimuliOhm', 1e6);

rowList = [rowList, {...
'TDS Output ', ...
'Ant 1 = Signal, Ant 2 = Signal, Ant 3 = Signal, Stimuli = 100kohm', ...
'ID43 = Mode 4 (cal mode 0), TDS_1 = V1_DC', ...
'ID44 = Mode 4 (cal mode 0), TDS_2 = V2_DC', ...
'ID45 = Mode 4 (cal mode 0), TDS_3 = V3_DC' }];
exp(end+1) = create_calib_table_record(43, 4, 1, [1]);
exp(end+1) = create_calib_table_record(44, 4, 2, [2]);
exp(end+1) = create_calib_table_record(45, 4, 3, [3]);
exp = set_field(exp, 3, 'antennaSignals', [1 1 1]);
exp = set_field(exp, 3, 'stimuliOhm', 1e5);

res = bicas.sact.parse_testlogbook(rowList);

if ~isequaln(exp, res)
    error('FAIL')
end

end



function Record = create_calib_table_record(testIdNbr, muxMode, outputChNbr, inputChNbr, varargin)
% Convenience function for assigning components in struct array.
% NOTE: Does not assign all values, only some.
%
% NOTE: Arguments are deliberately in the same order as the quoted testlogbook*.txt printouts.
% PROPOSAL: Make into nested function
% PROPOSAL: Abolish.
if length(varargin) == 0
    acGain = NaN;
elseif length(varargin) == 1
    acGain = varargin{1};
else
    error('BICAS:create_calib_table_record:Assertion', 'Illegal number of arguments.')
end

% 'mebTempCelsius', [], ...
% 'data',           []
Record = struct(...
    'antennaSignals', [], ...
    'stimuliOhm',     [], ...
    'testIdNbr',      testIdNbr, ...
    'muxMode',        muxMode, ...
    'outputChNbr',    outputChNbr, ...
    'inputChNbr',    inputChNbr, ...
    'acGain',         acGain);
end



function StructArray = set_field(StructArray, nLastIndices, fieldName, fieldValue)
% Utility function for assigning a specified field in the last N components in a 1D struct array.

for i = numel(StructArray) + [(-nLastIndices+1) : 0]
    StructArray(i).(fieldName) = fieldValue;
end
end



