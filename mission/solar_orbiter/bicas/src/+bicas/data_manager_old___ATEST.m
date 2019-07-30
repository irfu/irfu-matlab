function data_manager_old___ATEST
% data_manager_old_TST - Attempt at automated test code for data_manager_old.
%
% NOTE 2019-07-24: Does not work
% 
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-10-05
%
% NOTE: There can be NaN in individual samples.
%    Ex: Mysterious signal 7: ROC-SGSE_L2S_RPW-LFR-SURV-SWF-E_V01___OUTPUT_AV.cdf   
%       V: (35, 1477,2), (53, 500,3) (32, 1477,3)
%       E:  (32+3, 1477, 1)
%           (53+3,  500, 2)
%           (53+3,  500, 3)
%           (32+3, 1477, 2)


    % PROPOSAL: Create long complex test data from which one then deletes different subsets of "records", perhaps based on time.
    % QUESTION: How handle very similar datasets?
    % TODO: Do for mysterious signal i.e. V02, SURV-CWF, SURV-SWF!!
    
    test_1
end



function test_1
    % Testing V01_ROC-SGSE_L2R_RPW-LFR-SURV-CWF --> V01_ROC-SGSE_L2S_RPW-LFR-SURV-CWF-E
    % MUX mode 0
    %
    % NOTE: Does not test IBIAS.

    % PROPOSAL: Use dm_utils.* functions (==>This code does not test them), but create test code for those functions
    %           specifically instead.
    %   Ex: convert_N_to_1_SPR_ACQUISITION_TIME
    %   Ex: convert_N_to_1_SPR_Epoch
    % TODO: At least two rows - Tests reshaping snapshot/record-->sample/record better.
    
    clear   % Remove all global variables, i.e. CONSTANTS.    
    global CONSTANTS
    CONSTANTS = bicas.constants();
    ACQUISITION_TIME_EPOCH_UTC = [2000,01,01, 12,00,00, 000,000,000];

    Dm = bicas.data_manager_old();
    
    V_seq1   = rand(1,672);
    V_seq2   = rand(1,672);
    E_seq1_1  = rand(1,672);   % ELECTRICAL_1, record 1
    E_seq2_1  = rand(1,672);
    E_seq1_2  = rand(1,672);
    E_seq2_2  = rand(1,672);
    NaN_seq = ones(1,672) * NaN;
    
    HK = struct();
    HK.Epoch = int64(linspace(double(spdfparsett2000('2016-10-14 00:00:00')), double(spdfparsett2000('2016-10-14 01:00:00')), 2));
    HK.ACQUISITION_TIME    = uint32([0, 0; 1, 0]);
    HK.HK_BIA_MODE_MUX_SET = [0; 0];
    HK.HK_BIA_DIFF_GAIN    = [0; 0];
    
    SCI = struct();    
    SCI.Epoch = spdfparsett2000('2016-10-14 00:30:00')+int64([0;1e9]);
    SCI.ACQUISITION_TIME = uint32([0, 32767; 1, 0]);
    SCI.QUALITY_FLAG = [];
    SCI.QUALITY_BITMASK = [];
    SCI.POTENTIAL  = [V_seq1, V_seq2];
    %SCI.ELECTRICAL(:,:,1) = [E_seq1_1];
    %SCI.ELECTRICAL(:,:,2) = [E_seq2_1];
    SCI.ELECTRICAL(:,:,1) = [E_seq1_1; E_seq1_2];
    SCI.ELECTRICAL(:,:,2) = [E_seq2_1; E_seq2_2];
    SCI.FREQ = [2; 2];
    SCI.R0 = [0; 0];
    SCI.R1 = [0; 0];
    SCI.R2 = [0; 0];

    SCI_out = struct;
    SCI_out.Epoch            = bicas.dm_utils.convert_N_to_1_SPR_Epoch(          SCI.Epoch,            672, [256; 256]);
    SCI_out.ACQUISITION_TIME = bicas.dm_utils.convert_N_to_1_SPR_ACQUISITION_TIME(SCI.ACQUISITION_TIME, 672, [256; 256], ACQUISITION_TIME_EPOCH_UTC);
    
    SCI_out.V = [V_seq1'*17, NaN_seq', NaN_seq'; V_seq2'*17, NaN_seq', NaN_seq'];
    SCI_out.E = [NaN_seq',  NaN_seq', NaN_seq'];
    SCI_out.EAC = [E_seq1_1', E_seq1_1'+E_seq2_1', E_seq2_1'] / 5;

%     Dm.set_elementary_input_process_data('V01_ROC-SGSE_HK_RPW-BIA', HK);
    Dm.set_elementary_input_process_data('V02_ROC-SGSE_HK_RPW-BIA', HK);    
    Dm.set_elementary_input_process_data('V01_ROC-SGSE_L2R_RPW-LFR-SURV-CWF', SCI);
    %DM.set_elementary_input_process_data('V02_ROC-SGSE_L2R_RPW-LFR-SURV-CWF', SCI);
%     SCI_out_result = Dm.get_process_data_recursively('V01_ROC-SGSE_L2S_RPW-LFR-SURV-CWF-E', 'LFR-SURV-CWF-E_V01-V01');
    SCI_out_result = Dm.get_process_data_recursively('V03_ROC-SGSE_L2S_RPW-LFR-SURV-CWF-E');
    
    
    
    epsilon = 1e-13;
    fn_list = fieldnames(SCI_out)';
    for fn = fn_list
        A = SCI_out.(fn{1});
        B = SCI_out_result.(fn{1});
        %if bicas.utils.equals_tolerance(A,B, epsilon)
        if EJ_library.utils.equals_recursive(A,B, 'epsilon', epsilon)
            %fprintf('%s - Matches\n', fn{1})
        else
            fprintf('%s - NO MATCH\n', fn{1})
            A(1:3, :)
            B(1:3, :)
            keyboard
        end
    end
    
    %SCI_out.Epoch == SCI_out_result
    %SCI_out.Epoch(end-10:end)
    %SCI_out_result.Epoch(end-10:end)
    
end



