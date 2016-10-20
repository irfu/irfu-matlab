function data_manager_TEST
% data_manager_TST - Automated test code for data_manager.
% 
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-10-05
%

    % PROPOSAL: Create long complex test data from which one then deletes different subsets of "records", perhaps based on time.
    % QUESTION: How handle very similar datasets?
    % TODO: Do for mysterious signal i.e. V02, SURV-CWF, SURV-SWF!!
    
    test_1
end



%function test_2
%end



function test_1
    % Testing L2R_LFR-SURV-CWF_V01 --> L2S_LFR-SURV-CWF-E_V01
    % MUX mode 0
    %
    % NOTE: Does not test IBIAS.

    % PROPOSAL: Use dm_utils.* functions (==>This code does not test them), but create test code for those functions
    %           specifically instead.
    %   Ex: ACQUISITION_TIME___expand_to_sequences
    %   Ex: tt2000___expand_to_sequences
    % TODO: At least two rows - Tests reshaping snapshot/record-->sample/record better.
    
    clear   % Remove all global variables, i.e. CONSTANTS.    
    global CONSTANTS
    CONSTANTS = bicas.constants('');

    DM = bicas.data_manager();
    
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
    SCI_out.Epoch            = bicas.dm_utils.tt2000___expand_to_sequences(          SCI.Epoch,            672, [256; 256]);
    SCI_out.ACQUISITION_TIME = bicas.dm_utils.ACQUISITION_TIME___expand_to_sequences(SCI.ACQUISITION_TIME, 672, [256; 256]);
    
    SCI_out.V = [V_seq1'*17, NaN_seq', NaN_seq'; V_seq2'*17, NaN_seq', NaN_seq'];
    SCI_out.E = [NaN_seq',  NaN_seq', NaN_seq'];
    SCI_out.EAC = [E_seq1_1', E_seq1_1'+E_seq2_1', E_seq2_1'] / 5;

    DM.set_elementary_input_process_data('HK_BIA_V01', HK);
    DM.set_elementary_input_process_data('L2R_LFR-SURV-CWF_V01', SCI);
    %DM.set_elementary_input_process_data('L2R_LFR-SURV-CWF_V02', SCI);
    SCI_out_result = DM.get_process_data_recursively('L2S_LFR-SURV-CWF-E_V01', 'LFR-SURV-CWF-E_V01-V01');
    
    epsilon = 1e-13;
    fn_list = fieldnames(SCI_out)';
    for fn = fn_list
        A = SCI_out.(fn{1});
        B = SCI_out_result.(fn{1});
        if bicas.utils.equals_tolerance(A,B, epsilon)
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



