% Return MATLAB structure that "exactly" corresponds to the S/W descriptor
% specified by the RCS ICD. Return result can be used for generating a string
% that can be printed.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created ~2016-06-01
%
%
% RETURN VALUE
% ============
% JsonSwd : See bicas.utils.JSON_object_str for the exact format.
%
%
% IMPLEMENTATION NOTE
% ===================
% The values are derived entirely from BICAS constants. The structure is NOT
% directly incorporated in the BICAS constants/settings for extra flexibility.
% Reasons for NOT putting the S/W descriptor structure inside the settings
% class:
% (1) Some of the S/W descriptor variables have vague or misleading names
%     ("name", "dataset versions", "dataset IDs") which would (reasonably) have
%     to be represented by MATLAB variables with the same names.
% (2) Some of the S/W descriptor variables are grouped in a way which does not
%     fit the rest of the code (e.g. modes[i].outputs.(output_XX).release in the
%     S/W descriptor structure).
% (3) Some of the S/W descriptor values are really structure field NAMES, but
%     would be better as structure field VALUES (e.g. input CLI parameter,
%     output JSON identifier string).
% (4) The constants structure would become dependent on the format of the S/W
%     descriptor structure. The latter might change in the future, or be
%     misunderstood in the present (i.e. be changed).
% (5) Some S/W descriptor values repeat or can be derived from other constants
%     (e.g. author, contact, institute, output_XX.release.file).
% (6) It is easier to add automatic checks on the S/W descriptor in the code
%     that derives it.
%
function JsonSwd = get_SWD(SwmList)
%
% PROPOSAL: Have this function add the prefix "input_" to all modes[].inputs[].input (in SWD)
% and only store the "CLI_PARAMETER_suffix" in the BICAS constants structure instead.
%   CON: The main function also uses the full CLI parameter. It too would have to add the prefix.
%        ==> The same "input_" prefix is specified in two places ==> Bad practice.
%
% PROPOSAL: Validation of the information.
%   PROPOSAL: Separate S/W descriptor validation code which always runs when BICAS is run.
%       PRO: Implementation (where the values are taken from) may change.
%   PROPOSAL: Use the JSON schema in the docs to verify the S/W descriptor automatically.
%       NOTE: Requires there to be code for validating.
%   PROPOSAL: All versions, CLI parameters, dates, dataset IDs on the right format.
%   PROPOSAL: No mode name, dataset ID doubles.
%   PROPOSAL: Check that there are no additional structure fields.
%   PROPOSAL: Check that paths, filenames are valid.
%   QUESTION: How handle overlap validation of BICAS constants, since all values come from there?
%   NOTE: Already implicitly checks that all the needed fields exist (since they are read here).
%   NOTE: Checks on the main constants structure will (can) only happen if this file is executed, not if
%         the S/W as a whole is (by default).

assert(isvector(SwmList))
assert(isa(SwmList, 'bicas.swm.SoftwareMode'))

% Variable naming convention:
% ---------------------------
% Swd = S/W descriptor
%   The MATLAB structure which is used for producing the S/W descriptor
%   (SWD) JSON object string. Its fields (field names) should NOT follow
%   variable naming conventions since they determine the JSON object string
%   which must follow the RCS ICD.
SwdMetadataMap = bicas.const.SWD_METADATA;

JsonSwd = [];
JsonSwd.identification.project     = SwdMetadataMap('SWD.identification.project');
JsonSwd.identification.name        = SwdMetadataMap('SWD.identification.name');
JsonSwd.identification.identifier  = SwdMetadataMap('SWD.identification.identifier');
JsonSwd.identification.description = SwdMetadataMap('SWD.identification.description');
JsonSwd.identification.icd_version = SwdMetadataMap('SWD.identification.icd_version');

JsonSwd.release.version            = SwdMetadataMap('SWD.release.version');
JsonSwd.release.date               = SwdMetadataMap('SWD.release.date');
JsonSwd.release.author             = SwdMetadataMap('SWD.release.author');
JsonSwd.release.contact            = SwdMetadataMap('SWD.release.contact');
JsonSwd.release.institute          = SwdMetadataMap('SWD.release.institute');
JsonSwd.release.modification       = SwdMetadataMap('SWD.release.modification');
% RCS ICD 00037 iss1/rev2, draft 2019-07-11: Optional.
JsonSwd.release.source             = SwdMetadataMap('SWD.release.source');

JsonSwd.environment.executable     = SwdMetadataMap('SWD.environment.executable');
% RCS ICD 00037 iss1/rev2, draft 2019-07-11: Optional.
JsonSwd.environment.configuration  = bicas.const.DEFAULT_CONFIG_FILE_RELATIVE_PATH;

JsonSwd.modes = {};
for i = 1:length(SwmList)
  JsonSwd.modes{end+1} = generate_SWD_mode(SwmList(i));
end

end



function JsonSwdMode = generate_SWD_mode(Swm)
assert(isa(Swm, 'bicas.swm.SoftwareMode'))

JsonSwdMode.name    = Swm.cliOption;
JsonSwdMode.purpose = Swm.swdPurpose;

JsonSwdMode.inputs = [];
for i = 1:length(Swm.inputsList)
  InputDataset = Swm.inputsList(i);

  JsonInput = [];
  JsonInput.identifier = InputDataset.dsi;
  JsonSwdMode.inputs.(InputDataset.cliOptionHeaderBody) = JsonInput;
end

JsonSwdMode.outputs = {};
for i = 1:length(Swm.outputsList)
  OutputDataset = Swm.outputsList(i);

  JsonOutput = [];
  JsonOutput.identifier  = OutputDataset.dsi;
  JsonOutput.name        = OutputDataset.swdName;
  JsonOutput.description = OutputDataset.swdDescription;
  JsonOutput.level       = OutputDataset.datasetLevel;
  JsonOutput.template    = bicas.get_master_CDF_filename(...
    OutputDataset.dsi, ...
    OutputDataset.skeletonVersion);    % RCS ICD 00037 iss1/rev2, draft 2019-07-11: Optional.
  JsonSwdMode.outputs.(OutputDataset.cliOptionHeaderBody) = JsonOutput;
end

end
