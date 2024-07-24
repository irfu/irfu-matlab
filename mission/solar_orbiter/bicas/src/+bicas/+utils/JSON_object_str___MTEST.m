%
% Informal manual test code.
%
% NOTE: The example s/w descriptor is outdated, but that does not really matter
% for the purpose of manually testing the JSON format.
%
function JSON_object_str___MTEST
% PROPOSAL: Delete.
%   CON: Useful for verifying BICAS's old implementation (obsoleted but only
%        commented out).
%   CON-PROPOSAL: Delete when deleting commented-out old implementation.
%     NOTE: There is no automated test code.

obj = define_descriptor1();
str = bicas.utils.JSON_object_str(obj, 4);
fprintf('-------------------------------------------\n');
fprintf(str);



%==========================================================================
  function obj = define_descriptor3()
    obj = struct('output_cdf1', 'output_filename1.cdf', 'output_cdf2', 'output_filename2.cdf');
  end
%==========================================================================
  function obj = define_descriptor2()
    obj = {};
    obj{1} = struct('output_cdf1', 'output_filename1.cdf');
    obj{2} = struct('output_cdf2', 'output_filename2.cdf');
  end
%==========================================================================
  function obj = define_descriptor1()
    ERIK_P_G_JOHANSSON = 'Erik P G Johansson';

    obj = struct();

    identification.project     = 'ROC-SGSE';
    identification.name        = 'BICAS (temporary name)';
    identification.identifier  = 'ROC-SGSE-BICAS';    % Temporary
    identification.description = 'BIAS calibration software (temporary description)';
    obj.identification = identification;

    release.version = '0.0.1';
    release.date         = '2016-05-19';
    release.author       = ERIK_P_G_JOHANSSON;
    release.contact      = 'erik.johansson@irfu.se';
    release.institute    = 'IRF-U';
    release.modification = 'None (Initial release)';
    obj.release = release;

    environment.executable = 'bin/bicas';
    obj.environment = environment;

    obj.modes = [];

    mode = [];
    mode.name = 'testmode1';
    mode.purpose =  'Mode 1 (temporary purpose description)';
    mode.inputs.input_SCI.identifier = 'ROC-SGSE_L2R_RPW-LFR-SBM1-CWF';
    mode.inputs.input_SCI.version    = '01';
    mode.inputs.input_HK.identifier  = 'ROC-SGSE_HK_RPW-BIA';
    mode.inputs.input_HK.version     = '01';

    mode.outputs.output_SCI.identifier  = 'ROC-SGSE_L2S_RPW-BIA-xxxxx';
    mode.outputs.output_SCI.name        = 'xxxxx (temporary name)';
    mode.outputs.output_SCI.description = 'Contains xxxxx (temporary description)';
    mode.outputs.output_SCI.level       = 'L2S';
    mode.outputs.output_SCI.release.author       = ERIK_P_G_JOHANSSON;
    mode.outputs.output_SCI.release.date         = '2016-05-19';
    mode.outputs.output_SCI.release.version      = '01';
    mode.outputs.output_SCI.release.contact      = 'erik.johansson@irfu.se';
    mode.outputs.output_SCI.release.institute    = 'IRF-U';
    mode.outputs.output_SCI.release.modification = 'None (initial release)';
    obj.modes{end+1} = mode;

    mode = [];
    mode.name = 'testmode2';
    mode.purpose =  'Mode 2 (temporary purpose description)';
    mode.inputs.input_SCI.identifier = 'ROC-SGSE_L2R_RPW-TDS-SURV-RSWF';
    mode.inputs.input_SCI.version    = '01';
    mode.outputs.output_SCI.identifier  = 'ROC-SGSE_L2S_RPW-BIA-xxxxx';
    mode.outputs.output_SCI.name        = 'xxxxx (temporary name)';
    mode.outputs.output_SCI.description = 'Contains xxxxx (temporary description)';
    mode.outputs.output_SCI.level       = 'L2S';
    mode.outputs.output_SCI.release.author       = ERIK_P_G_JOHANSSON;
    mode.outputs.output_SCI.release.date         = '2016-05-19';
    mode.outputs.output_SCI.release.version      = '01';
    mode.outputs.output_SCI.release.contact      = 'erik.johansson@irfu.se';
    mode.outputs.output_SCI.release.institute    = 'IRF-U';
    mode.outputs.output_SCI.release.modification = 'None (initial release)';
    obj.modes{end+1} = mode;
  end
end
