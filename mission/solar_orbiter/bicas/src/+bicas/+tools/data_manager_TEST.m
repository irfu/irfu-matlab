% data_manager_TST - Automated test code for data_manager.
% 
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-10-05
%
function data_manager_TEST

    DM = bicas.data_manager();
    
    HK = struct();
    SCI = struct();
    
    DM.set_elementary_input_process_data('HK_BIA_V01', HK);
    DM.set_elementary_input_process_data('L2R_LFR-SURV-CWF_V01', SCI);
    
    %DM.get_process_data_recursively('L2S_LFR-SURV-CWF-E_V01', 'LFR-SURV-CWF-E_V01-V01')
end