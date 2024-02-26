%
% Test performance of
% bicas.tools.batch.autocreate_input_BPCIs().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function autocreate_input_BPCIs2___speedTest()

    % profile clear; profile on;

    nArray           = [];
    wallTimeSecArray = [];

    N_CALLS = 5;

    for n = round(logspace(log10(20), log10(1000), N_CALLS))
        % NOTE: test() prints log messages.
        wallTimeSec = test(n);

        % Log
        n
        wallTimeSec

        nArray(end+1)           = n;
        wallTimeSecArray(end+1) = wallTimeSec;
    end

    close all
    figure('WindowState', 'maximized')
    plot(nArray, wallTimeSecArray, 'o-')
    grid on

    % profile off; profile viewer
end



function wallTimeSec = test(n)
    DSI_1 = 'SOLO_L1R_RPW-LFR-SURV-CWF-E';
    DSI_2 = 'SOLO_L2_RPW-LFR-SURV-CWF-E';

    SWMP = bicas.tools.batch.TestSwmProcessing();

    SWM_1 = bicas.swm.SoftwareMode(...
        SWMP, 'CLI_SWM_1', 'SWD purpose', ...
        bicas.swm.InputDataset( 'cli_in',  DSI_1, 'IN_cdf'), ...
        bicas.swm.OutputDataset('cli_out', DSI_2, 'OUT_cdf', 'SWD ', 'SWD ', '02') ...
    );

    path1Ca = get_dataset_paths(n, 'solo_L1R_rpw-lfr-surv-cwf-e_%4i%02i%02i_V02.cdf');
    path2Ca = get_dataset_paths(n, 'solo_L2_rpw-lfr-surv-cwf-e_%4i%02i%02i_V01.cdf');

    InputDsmdArray               = solo.adm.paths_to_DSMD_array(path1Ca);
    preexistingOutputFilenamesCa = path2Ca;   % Collision with every output file.

    OUTPUT_DIR    = '/tmp/nonexisting';    % Should be irrelevant.
    OUTPUT_ISCDAG = false;                 % Must match preexistingOutputFilenamesCa.
    CURRENT_DATASET_EXTENSION_DAYS = 0;

    PreexistingOutputLvDsmdArray = solo.adm.paths_to_DSMD_array(preexistingOutputFilenamesCa);

    get_BPCI_output_path_fh = ...
        @(outputDsi, BpciInputDsmdArray) ( ...
            bicas.tools.batch.get_BPCI_output_path2( ...
                BpciInputDsmdArray, PreexistingOutputLvDsmdArray, ...
                outputDsi, ...
                'HIGHEST_USED', ...
                OUTPUT_DIR, OUTPUT_ISCDAG ...
            ) ...
        );

    % ===============
    % RUN TESTED CODE
    % ===============
    t = tic();
    BpciInputArray = bicas.tools.batch.autocreate_input_BPCIs(...
        InputDsmdArray, get_BPCI_output_path_fh, [SWM_1], ...
        CURRENT_DATASET_EXTENSION_DAYS);
    wallTimeSec = toc(t);

    assert(numel(BpciInputArray) == numel(InputDsmdArray))

end



% Generate n dataset paths.
%
% NOTE: It is important that code generates correct dataset filenames which can be
% recognized as datasets. Otherwise the tested code will ignore them, and
% presumable run faster.
function filenameCa = get_dataset_paths(n, patternStr)
    dt = datetime('2020-01-01') + caldays([1:n]-1)';

    filenameCa = cell(n, 1);
    for i = 1:n
        assert(dt(i).Year < 10000)   % Year must be four digits.

        filenameCa{i} = sprintf(patternStr, dt(i).Year, dt(i).Month, dt(i).Day);
    end
end
