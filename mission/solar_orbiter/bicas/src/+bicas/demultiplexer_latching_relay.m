%
% Return the value for the BIAS demultiplexer latching relay which, together with the demultiplexer mode, determines how
% BLTS channels and ASR signals are related.
%
% The latching relay state should ideally be obtained from telecommanding just like the bias currents, but this has not
% been implemented into the ROC pipeline and is not planned to be since the latching relay state is not expected to be
% changed unless a hardware failure (probe failure) occurs (2019-11-19). Therefore, the latching relay state is
% hardcoded until further notice
%
%
% ARGUMENTS AND RETURN VALUE
% ==========================
% Epoch 
% dlrUsing12 : 0/1, true/false. Array same size as Epoch.
%               False=0 = Using diffs V13_DC, V13_AC
%               True =1 = Using diffs V12_DC, V12_AC
%
%
% DEFINITIONS
% ===========
% DLR : Demultiplexer Latching Relay
% See bicas.calib.
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2019-11-18
%
function dlrUsing12 = demultiplexer_latching_relay(Epoch)
% PROPOSAL: SETTING for overriding.

    bicas.proc_utils.assert_Epoch(Epoch)
    
    dlrUsing12 = ones(size(Epoch));
    
end
