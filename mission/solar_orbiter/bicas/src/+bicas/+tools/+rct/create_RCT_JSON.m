%
% Create RCT JSON file, associated with BICAS RCTs.
%
% NOTE: The RCT JSON filename should be static and is therefore hard-coded.
% NOTE: Will overwrite old RCT JSON file.
%
%
% SPECIFICATION OF JSON FILE
% ==========================
% """"""""
% Please find attached an "dummy" example of what the RCT catalog JSON file
% could look like.
%
% Main specifications in mind are:
%
% 	- Format is JSON
% 	- One file per team
% 	- Static name. I propose something simple —> lfr/bia/scm/tds/thr_rct_validity.json
% 	- File must be delivered by teams with any new RCT file in the cal/
%     subfolder in the LESIA sftp server
%
%     And to make it simple as possible:
%
%     - No header
%     - One field (JSON object) per RCT file (e.g.
%       "SOLO_CAL_RCT-LFR-SCM_V20190123171020.cdf")
%     - Each RCT field gives at least one validity time range (i.e.,
%       "validity_start" and "validity_end").
%
%
% Let me know if you have comments or suggestions of improvement (if not I will
% assume that you approve the proposed content).
% """""""" /Xavier Bonnin e-mail 2020-06-12
%
% lfr_rct_validity.json:
% {
%     "SOLO_CAL_RCT-LFR-VHF_V20190123171020.cdf":[
%             {
%                 "validity_start":"2020-02-10T00:00:00Z",
%                 "validity_end": "2020-05-01T00:00:00Z"
%             },
%             {
%                 "validity_start": "2021-05-01T00:00:00Z",
%                 "validity_end": "9999-01-01T00:00:00Z"
%             }
%     ],
%     "SOLO_CAL_RCT-LFR-SCM_V20190123171020.cdf":[
%             {
%                 "validity_start": "2020-02-10T00:00:00Z",
%                 "validity_end": "9999-01-01T00:00:00Z"
%             }
%         ],
%     "SOLO_CAL_RCT-LFR-BIAS_V20190123171020.cdf":[
%             {
%                 "validity_start": "2020-02-10T00:00:00Z",
%                 "validity_end": "9999-01-01T00:00:00Z"
%             }
%         ]
% }
% /Xavier Bonnin e-mail 2020-06-12
%
%
% RETURN VALUE
% ============
% rctJsonPath : Path to the file created. Useful for printing log messages.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-06-24.
%
function rctJsonPath = create_RCT_JSON(destDir, biasRctFilename, beginDt, endDt)
% PROPOSAL: Do not use bicas.utils.JSON_object_str(). Use MATLAB's own support
%           for JSON files: jsonencode().
DT_FORMAT_STR = 'yyyy-MM-dd''T''HH:mm:ss''Z''';

beginStr = char(datetime(beginDt, 'Format', DT_FORMAT_STR));
endStr   = char(datetime(endDt,   'Format', DT_FORMAT_STR));


RCT_JSON_FILENAME = 'bias_rct_validity.json';

% NOTE: Cell array of struct, to conform with XB's format above.
JsonObj = containers.Map();
JsonObj(biasRctFilename) = { ...
  struct(...
    'validity_start', beginStr, ...
    'validity_end',   endStr ...
  ) ...
};

str = bicas.utils.JSON_object_str(JsonObj, 4);
%fprintf(str);    % DEBUG

rctJsonPath = fullfile(destDir, RCT_JSON_FILENAME);
irf.fs.write_file(rctJsonPath, uint8(str(:)));

end