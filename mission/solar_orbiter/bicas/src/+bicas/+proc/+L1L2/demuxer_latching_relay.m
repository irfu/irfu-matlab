%
% Return the value for the BIAS demultiplexer latching relay which, together
% with the demultiplexer mode, determines how BLTS channels and ASR signals are
% related.
%
% The latching relay state should ideally be obtained from telecommanding just
% like the bias currents, but this has not been implemented into the ROC
% pipeline and is not planned to be since the latching relay state is not
% expected to be changed unless a hardware failure (probe failure) occurs
% (2019-11-19). Therefore, the latching relay state is hard-coded until further
% notice.
%
% NOTE: BIAS HK contains the value in HK_BIA_MODE_DIFF_PROBE (presumably). See
% BIAS specification, section "3.4.4.14 MODE", "Data D3 = Diff probe 1&2(0),
% Diff probe 1&3(1)"
%
% 2023-07-28: YK has asked ROC/Diane Berard who has asked SOC to change the DLR
% value. "will likely take several weeks".
% Diane Berard e-mail 2023-08-29: The latching relay did change on "August 21st
% at 12:03:04 in flight". HK_BIA_MODE_DIFF_PROBE should indeed (very likely) be
% the latching relay.
%
%
% ARGUMENTS AND RETURN VALUE
% ==========================
% Epoch 
% dlrFpa
%       Logical. FPA same size as Epoch.
%       NOTE: The interpretation of values follows the same convention as
%       BIAS HK ZV "HK_BIA_MODE_DIFF_PROBE".
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2019-11-18
%
function dlrFpa = demuxer_latching_relay(Epoch)
    bicas.utils.assert_ZV_Epoch(Epoch)

    dlrFpa = bicas.utils.FillPositionsArray(...
        false(size(Epoch)), 'fill positions', false(size(Epoch)));
end
