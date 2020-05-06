function res = get_data(varStr, Tint, mmsId)
%MMS.GET_DATA  Load a variable
%
%  res = MMS.GET_DATA(varStr, Tint, [mmsId])
%
%  mmsId=0 (DEFAULT) means 1:4
%
%  varStr is one of:
%  EPHEMERIS:
%      'R_gse', 'R_gsm', 'V_gse', 'V_gsm', ...
%      'tetra_quality'.
%  EDP:
%      'Phase_edp_fast_l2a', 'Phase_edp_slow_l2a',...
%      'Es12_dsl_edp_fast_l2a', 'Es12_dsl_edp_slow_l2a',...
%      'Es34_dsl_edp_fast_l2a', 'Es34_dsl_edp_slow_l2a',...
%      'Sdev12_edp_fast_l2a', 'Sdev12_edp_slow_l2a',...
%      'Sdev34_edp_fast_l2a', 'Sdev34_edp_slow_l2a',...
%      'Adcoff_edp_fast_l2a', 'Adcoff_edp_slow_l2a',...
%      'E_dsl_edp_brst_l2', 'E_dsl_edp_fast_l2',...
%      'E_gse_edp_brst_l2', 'E_gse_edp_fast_l2',...
%      'E2d_dsl_edp_brst_l2pre', 'E2d_dsl_edp_fast_l2pre',...
%      'E2d_dsl_edp_l2pre', 'E_dsl_edp_l2pre',...
%      'E_ssc_edp_brst_l1b', 'E_ssc_edp_fast_l1b', 'E_ssc_edp_slow_l1b'.
%  FPI IONS:
%     'Ni_fpi_brst_l2' (alias:'Ni_fpi_brst'), 'partNi_fpi_brst_l2', ...
%     'Ni_fpi_fast_l2', 'Ni_fpi_sitl', 'Ni_fpi_ql',...
%     'Ni_fpi_brst_l1b', 'Ni_fpi_fast_l1b', 'partNi_fpi_fast_l2',...
%     'Vi_dbcs_fpi_brst_l2' (alias:'Vi_dbcs_fpi_brst'), ...
%     'partVi_dbcs_fpi_brst_l2', 'Vi_dbcs_fpi_fast_l2',...
%     'Vi_gse_fpi_sitl', 'Vi_gse_fpi_ql', 'Vi_dbcs_fpi_fast_ql',...
%     'Vi_gse_fpi_brst_l1b', 'Vi_gse_fpi_fast_l1b',...
%     'Vi_gse_fpi_brst_l2', 'partVi_gse_fpi_brst_l2', ...
%     'Tsi_fpi_brst_l2' (alias:'Tsi_fpi_brst'), 'Tsi_fpi_fast_l2',... %scalar temperature
%     'Tperpi_fpi_brst_l2', 'Tparai_fpi_brst_l2',...
%     'Tsi_fpi_sitl', 'Tsi_fpi_ql',...
%     'Tsi_fpi_brst_l1b', 'Tsi_fpi_fast_l1b',...
%     'Ti_dbcs_fpi_brst_l2' (alias:'Ti_dbcs_fpi_brst'), ...
%     'partTi_dbcs_fpi_brst_l2', 'Ti_dbcs_fpi_fast_l2',...
%     'Ti_gse_fpi_sitl', 'Ti_gse_fpi_ql', 'Ti_dbcs_fpi_ql',...
%     'Ti_gse_fpi_brst_l1b', 'Ti_gse_fpi_fast_l1b',...
%     'Ti_gse_fpi_brst_l2', 'partTi_gse_fpi_brst_l2', ...
%     'Pi_dbcs_fpi_brst_l2' (alias:'Pi_dbcs_fpi_brst'), ...
%     'partPi_dbcs_fpi_brst_l2', 'Pi_dbcs_fpi_fast_l2',...
%     'Pi_gse_fpi_sitl' (alias:'Pi_gse_fpi_ql'),...
%     'Pi_gse_fpi_brst_l1b', 'Pi_gse_fpi_fast_l1b',...
%     'Pi_gse_fpi_brst_l2', 'partPi_gse_fpi_brst_l2', ...
%     'Enfluxi_fpi_fast_ql', 'Energyi_fpi_fast_ql',...
%     Skymaps:
%     'PDi_fpi_brst_l2', 'PDi_fpi_fast_l2'.
%  FPI ELECTRONS:
%     'Ne_fpi_brst_l2' (alias:'Ne_fpi_brst'), 'partNe_fpi_brst_l2', ...
%     'Ne_fpi_fast_l2', 'Ne_fpi_sitl', 'Ne_fpi_ql', 'partNe_fpi_fast_l2',...
%     'Ne_fpi_brst_l1b', 'Ne_fpi_fast_l1b',...
%     'Ve_dbcs_fpi_brst_l2' (alias:'Ve_dbcs_fpi_brst'), ...
%     'partVe_dbcs_fpi_brst_l2', 'Ve_dbcs_fpi_fast_l2',...
%     'Ve_gse_fpi_sitl', 'Ve_gse_fpi_ql',...
%     'Ve_gse_fpi_brst_l1b', 'Ve_gse_fpi_fast_l1b',...
%     'Ve_gse_fpi_brst_l2', 'partVe_gse_fpi_brst_l2', ...
%     'Tse_fpi_brst_l2' (alias:'Tse_fpi_brst'), 'Tse_fpi_fast_l2',... %scalar temperature
%     'Tperpe_fpi_brst_l2', 'Tparae_fpi_brst_l2',...
%     'Tse_fpi_sitl', 'Tse_fpi_ql',...
%     'Tse_fpi_brst_l1b', 'Tse_fpi_fast_l1b',...
%     Loads into tensor of order 2:
%     'Te_dbcs_fpi_brst_l2' (alias:'Te_dbcs_fpi_brst'), ...
%     'partTe_dbcs_fpi_brst_l2', 'Te_dbcs_fpi_fast_l2',...
%     'Te_gse_fpi_sitl', 'Te_gse_fpi_ql', 'Te_dbcs_fpi_ql',...
%     'Te_gse_fpi_brst_l1b', 'Te_gse_fpi_fast_l1b',...
%     'Te_gse_fpi_brst_l2', 'partTe_gse_fpi_brst_l2', ...
%     'Pe_dbcs_fpi_brst_l2' (alias:'Pe_dbcs_fpi_brst'), ...
%     'partPe_dbcs_fpi_brst_l2', 'Pe_dbcs_fpi_fast_l2',...
%     'Pe_gse_fpi_sitl' (alias:'Pe_gse_fpi_ql'),...
%     'Pe_gse_fpi_brst_l1b', 'Pe_gse_fpi_fast_l1b',...
%     'Pe_gse_fpi_brst_l2', 'partPe_gse_fpi_brst_l2', ...
%     'Enfluxe_fpi_fast_ql', 'Energye_fpi_fast_ql', ... 
%     Skymaps:
%     'PDe_fpi_brst_l2', 'PDe_fpi_fast_l2'.
%  FGM:
%     'B_gsm_fgm_srvy_l2' (aliases:'B_gsm_srvy_l2', 'B_gsm_srvy'),...
%     'B_gsm_fgm_brst_l2' (aliases:'B_gsm_brst_l2', 'B_gsm_brst'),...
%     'B_gse_fgm_srvy_l2' (aliases:'B_gse_srvy_l2', 'B_gse_srvy'),...
%     'B_gse_fgm_brst_l2' (aliases:'B_gse_brst_l2', 'B_gse_brst'),...
%     'B_bcs_fgm_srvy_l2' (aliases:'B_bcs_srvy_l2', 'B_bcs_srvy')...
%     'B_bcs_fgm_brst_l2' (aliases:'B_bcs_brst_l2', 'B_bcs_brst'),...
%     'B_dmpa_fgm_srvy_l2' (aliases:'B_dmpa_srvy_l2', 'B_dmpa_srvy'),...
%     'B_dmpa_fgm_brst_l2' (aliases:'B_dmpa_brst_l2', 'B_dmpa_brst').
%  FSM:
%     'B_gse_fsm_brst_l3'.
%  DFG & AFG
%     'B_gsm_dfg_srvy_l2pre',  'B_gsm_afg_srvy_l2pre',...
%     'B_gse_dfg_srvy_l2pre',  'B_gse_afg_srvy_l2pre',...
%     'B_dmpa_dfg_srvy_l2pre', 'B_dmpa_afg_srvy_l2pre',...
%     'B_bcs_dfg_srvy_l2pre',  'B_bcs_afg_srvy_l2pre',...
%     'B_dmpa_dfg_srvy_ql',    'B_dmpa_afg_srvy_ql'.
%  HPCA:
%     'Nhplus_hpca_srvy_l2', 'Nheplus_hpca_srvy_l2', 'Nheplusplus_hpca_srvy_l2', 'Noplus_hpca_srvy_l2',...
%     'Tshplus_hpca_srvy_l2', 'Tsheplus_hpca_srvy_l2', 'Tsheplusplus_hpca_srvy_l2', 'Tsoplus_hpca_srvy_l2',...
%     'Vhplus_dbcs_hpca_srvy_l2', 'Vheplus_dbcs_hpca_srvy_l2', 'Vheplusplus_dbcs_hpca_srvy_l2', 'Voplus_dbcs_hpca_srvy_l2',...
%     'Phplus_dbcs_hpca_srvy_l2', 'Pheplus_dbcs_hpca_srvy_l2', 'Pheplusplus_dbcs_hpca_srvy_l2', 'Poplus_dbcs_hpca_srvy_l2',...
%     'Thplus_dbcs_hpca_srvy_l2', 'Theplus_dbcs_hpca_srvy_l2', 'Theplusplus_dbcs_hpca_srvy_l2', 'Toplus_dbcs_hpca_srvy_l2',...
%     'Vhplus_gsm_hpca_srvy_l2', 'Vheplus_gsm_hpca_srvy_l2', 'Vheplusplus_gsm_hpca_srvy_l2', 'Voplus_gsm_hpca_srvy_l2',...
%     'Phplus_gsm_hpca_srvy_l2', 'Pheplus_gsm_hpca_srvy_l2', 'Pheplusplus_gsm_hpca_srvy_l2', 'Poplus_gsm_hpca_srvy_l2',...
%     'Thplus_gsm_hpca_srvy_l2', 'Theplus_gsm_hpca_srvy_l2', 'Theplusplus_gsm_hpca_srvy_l2', 'Toplus_gsm_hpca_srvy_l2',...
%     'Nhplus_hpca_sitl',
%     'Omnifluxoplus_hpca_brst_l2','Omnifluxhplus_hpca_brst_l2','Omnifluxheplus_hpca_brst_l2','Omnifluxheplusplus_hpca_brst_l2'
%  EDI:
%     'Flux-amb-pm2_edi_brst_l2'
%  ASPOC:
%     'aspoc_status'.
%
% Example:
%   Tint = irf.tint('2015-09-21T00:00:00Z/2015-09-21T17:00:00Z');
%   V1 = mms.get_data('V_gse',Tint,1); % SC GSE velocity for MMS1
%   R  = mms.get_data('R_gse',Tint);   % SC GSE position for all MMS SC

res = [];

if nargin<3, mmsId = 0; end
if isempty(intersect(mmsId,0:4))
  errS = ['invalid MMS ID: ' mmsId];
  irf.log('critical',errS); error(errS)
end
mmsIdS = num2str(mmsId);

if ~isa(Tint,'GenericTimeArray')
  errS = 'TINT must be of GenericTimeArray type';
  irf.log('critical',errS); error(errS)
elseif Tint.stop-Tint.start<=0
  errS = 'TINT duration is zero or negative';
  irf.log('critical',errS); error(errS)
end

vars = {'R_gse','R_gsm','V_gse','V_gsm',...
  'B_gsm_fgm_srvy_l2','B_gsm_srvy_l2','B_gsm_srvy','B_gsm_fgm_brst_l2','B_gsm_brst_l2','B_gsm_brst',...
  'B_gse_fgm_srvy_l2','B_gse_srvy_l2','B_gse_srvy','B_gse_fgm_brst_l2','B_gse_brst_l2','B_gse_brst',...
  'B_bcs_fgm_srvy_l2','B_bcs_srvy_l2','B_bcs_srvy','B_bcs_fgm_brst_l2','B_bcs_brst_l2','B_bcs_brst',...
  'B_dmpa_fgm_srvy_l2','B_dmpa_srvy_l2','B_dmpa_srvy','B_dmpa_fgm_brst_l2','B_dmpa_brst_l2','B_dmpa_brst',...
  'B_gsm_dfg_srvy_l2pre',...
  'B_gse_dfg_srvy_l2pre',...
  'B_dmpa_dfg_srvy_l2pre',...
  'B_bcs_dfg_srvy_l2pre',...
  'B_gsm_afg_srvy_l2pre',...
  'B_gse_afg_srvy_l2pre',...
  'B_dmpa_afg_srvy_l2pre',...
  'B_bcs_afg_srvy_l2pre',...
  'B_dmpa_dfg_srvy_ql','B_dmpa_afg_srvy_ql',...
  'dfg_ql_srvy','afg_ql_srvy',...
  'B_gse_scm_brst_l2',...
  'B_gse_fsm_brst_l3', ...
  'tetra_quality',...
  'Phase_edp_fast_l2a','Phase_edp_slow_l2a',...
  'Es12_dsl_edp_slow_l2a','Es34_dsl_edp_slow_l2a',...
  'Es12_dsl_edp_fast_l2a','Es34_dsl_edp_fast_l2a',...
  'Sdev12_edp_slow_l2a','Sdev34_edp_slow_l2a',...
  'Sdev12_edp_fast_l2a','Sdev34_edp_fast_l2a',...
  'Adcoff_edp_fast_l2a','Adcoff_edp_slow_l2a',...
  'E_dsl_edp_brst_l2','E_dsl_edp_fast_l2','E_dsl_edp_brst_ql','E_dsl_edp_fast_ql','E_dsl_edp_slow_l2',...
  'E_gse_edp_brst_l2','E_gse_edp_fast_l2','E_gse_edp_brst_ql','E_gse_edp_fast_ql','E_gse_edp_slow_l2',...
  'E2d_dsl_edp_brst_l2pre','E2d_dsl_edp_fast_l2pre','E2d_dsl_edp_brst_ql','E2d_dsl_edp_fast_ql',...
  'E2d_dsl_edp_l2pre','E2d_dsl_edp_fast_l2pre','E2d_dsl_edp_brst_l2pre','E2d_dsl_edp_slow_l2pre',...
  'E_dsl_edp_l2pre','E_dsl_edp_fast_l2pre','E_dsl_edp_brst_l2pre','E_dsl_edp_slow_l2pre',...
  'E_ssc_edp_brst_l2a','E_ssc_edp_fast_l2a','E_ssc_edp_slow_l2a',...
  'E_ssc_edp_brst_l1b','E_ssc_edp_fast_l1b','E_ssc_edp_slow_l1b',...
  'Epar_edp_l2','Epar_edp_brst_l2','Epar_edp_fast_l2',...
  'V_edp_brst_l1b','V_edp_fast_l1b','V_edp_slow_l1b','V_edp_fast_sitl','V_edp_slow_sitl'...
  'V_edp_fast_l2','V_edp_brst_l2','V_edp_slow_l2',...
  'V6_edp_fast_l2','V6_edp_brst_l2','V6_edp_slow_l2',...
  'Vpsp_edp_fast_l2','Vpsp_edp_brst_l2','Vpsp_edp_slow_l2',...
  'Vi_dbcs_fpi_brst_l2', 'Vi_dbcs_fpi_brst', 'Vi_dbcs_fpi_fast_l2',...
  'Vi_gse_fpi_sitl', 'Vi_gse_fpi_ql', 'Vi_dbcs_fpi_ql', 'Vi_gse_fpi_fast_l2', ...
  'Vi_gse_fpi_brst_l1b','Vi_gse_fpi_brst_l2', 'partVi_gse_fpi_brst_l2', 'Vi_gse_fpi_fast_l1b',...
  'Ve_dbcs_fpi_brst_l2','Ve_dbcs_fpi_brst', 'Ve_dbcs_fpi_ql', 'Ve_dbcs_fpi_fast_l2',...
  'Ve_gse_fpi_sitl', 'Ve_gse_fpi_ql','Ve_gse_fpi_fast_l2',...
  'Ve_gse_fpi_brst_l1b','Ve_gse_fpi_fast_l1b', ...
  'Ve_gse_fpi_brst_l2', 'partVe_gse_fpi_brst_l2', ...
  'Ni_fpi_brst_l2', 'partNi_fpi_brst_l2', ...
  'Ni_fpi_brst','Ni_fpi_fast_l2','partNi_fpi_fast_l2',...
  'Ni_fpi_sitl','Ni_fpi_ql',...
  'Ni_fpi_brst_l1b','Ni_fpi_fast_l1b',...
  'Enfluxi_fpi_fast_ql', 'Enfluxe_fpi_fast_ql', ...
  'Energyi_fpi_fast_ql', 'Energye_fpi_fast_ql', ...
  'Ne_fpi_brst_l2', 'partNe_fpi_brst_l2', ...
  'Ne_fpi_brst','Ne_fpi_fast_l2','partNe_fpi_fast_l2',...
  'Ne_fpi_sitl','Ne_fpi_ql',...
  'Ne_fpi_brst_l1b','Ne_fpi_fast_l1b',...
  'Pe_fpi_ql','Pe_fpi_brst','Pe_fpi_brst_l2','Pi_fpi_brst_l2',...
  'Tsi_fpi_brst_l2','Tsi_fpi_brst','Tsi_fpi_fast_l2',...
  'Tperpi_fpi_brst_l2', 'Tparai_fpi_brst_l2', ...
  'partTperpi_fpi_brst_l2', 'partTparai_fpi_brst_l2', ...
  'Tsi_fpi_sitl','Tsi_fpi_ql',...
  'Tsi_fpi_brst_l1b','Tsi_fpi_fast_l1b',...
  'Tse_fpi_brst_l2','Tse_fpi_brst','Tse_fpi_fast_l2',...
  'Tperpe_fpi_brst_l2', 'Tparae_fpi_brst_l2', ...
  'partTperpe_fpi_brst_l2', 'partTparae_fpi_brst_l2', ...
  'Tse_fpi_sitl','Tse_fpi_ql',...
  'Tse_fpi_brst_l1b','Tse_fpi_fast_l1b',...
  'Ti_dbcs_fpi_brst_l2', 'partTi_dbcs_fpi_brst_l2', ...
  'Ti_dbcs_fpi_brst','Ti_dbcs_fpi_fast_l2',...
  'Ti_gse_fpi_sitl','Ti_gse_fpi_ql','Ti_dbcs_fpi_ql',...
  'Ti_gse_fpi_brst_l1b','Ti_gse_fpi_fast_l1b',...
  'Ti_gse_fpi_brst_l2', 'partTi_gse_fpi_brst_l2', ...
  'Te_dbcs_fpi_brst_l2', 'partTe_dbcs_fpi_brst_l2', ...
  'Te_dbcs_fpi_brst','Te_dbcs_fpi_fast_l2',...
  'Te_gse_fpi_sitl','Te_gse_fpi_ql','Te_dbcs_fpi_ql',...
  'Te_gse_fpi_brst_l1b','Te_gse_fpi_fast_l1b',...
  'Te_gse_fpi_brst_l2', 'partTe_gse_fpi_brst_l2', ...
  'Pi_dbcs_fpi_brst_l2', 'partPi_dbcs_fpi_brst_l2', ...
  'Pi_dbcs_fpi_brst','Pi_dbcs_fpi_fast_l2',...
  'Pi_gse_fpi_sitl','Pi_gse_fpi_ql',...
  'Pi_gse_fpi_brst_l1b','Pi_gse_fpi_fast_l1b',...
  'Pi_gse_fpi_brst_l2', 'partPi_gse_fpi_brst_l2',...
  'Pe_dbcs_fpi_brst_l2', 'partPe_dbcs_fpi_brst_l2', ...
  'Pe_dbcs_fpi_brst', 'Pe_dbcs_fpi_fast_l2',...
  'Pe_gse_fpi_sitl','Pe_gse_fpi_ql',...
  'Pe_gse_fpi_brst_l1b','Pe_gse_fpi_fast_l1b',...
  'Pe_gse_fpi_brst_l2', 'partPe_gse_fpi_brst_l2', ...
  'PDe_fpi_brst_l2','PDi_fpi_brst_l2',...
  'PDERRe_fpi_brst_l2','PDERRi_fpi_brst_l2',...
  'PDe_fpi_fast_l2','PDi_fpi_fast_l2',...
  'PDERRe_fpi_fast_l2','PDERRi_fpi_fast_l2',...
  'Flux-amb-pm2_edi_brst_l2','Flux-err-amb-pm2_edi_brst_l2',...
  'E_dsl_edi_fast_l2','E_dsl_edi_srvy_l2',...
  'Nhplus_hpca_srvy_l2','Nheplus_hpca_srvy_l2','Nheplusplus_hpca_srvy_l2','Noplus_hpca_srvy_l2',...
  'Tshplus_hpca_srvy_l2','Tsheplus_hpca_srvy_l2','Tsheplusplus_hpca_srvy_l2','Tsoplus_hpca_srvy_l2',...
  'Vhplus_dbcs_hpca_srvy_l2','Vheplus_dbcs_hpca_srvy_l2','Vheplusplus_dbcs_hpca_srvy_l2','Voplus_dbcs_hpca_srvy_l2',...
  'Phplus_dbcs_hpca_srvy_l2','Pheplus_dbcs_hpca_srvy_l2','Pheplusplus_dbcs_hpca_srvy_l2','Poplus_dbcs_hpca_srvy_l2',...
  'Thplus_dbcs_hpca_srvy_l2','Theplus_dbcs_hpca_srvy_l2','Theplusplus_dbcs_hpca_srvy_l2','Toplus_dbcs_hpca_srvy_l2',...
  'Vhplus_gsm_hpca_srvy_l2','Vheplus_gsm_hpca_srvy_l2','Vheplusplus_gsm_hpca_srvy_l2','Voplus_gsm_hpca_srvy_l2',...
  'Phplus_gsm_hpca_srvy_l2','Pheplus_gsm_hpca_srvy_l2','Pheplusplus_gsm_hpca_srvy_l2','Poplus_gsm_hpca_srvy_l2',...
  'Thplus_gsm_hpca_srvy_l2','Theplus_gsm_hpca_srvy_l2','Theplusplus_gsm_hpca_srvy_l2','Toplus_gsm_hpca_srvy_l2',...
  'Nhplus_hpca_srvy_l1b','Nheplus_hpca_srvy_l1b','Nheplusplus_hpca_srvy_l1b','Noplus_hpca_srvy_l1b',...
  'Tshplus_hpca_srvy_l1b','Tsheplus_hpca_srvy_l1b','Tsheplusplus_hpca_srvy_l1b','Tsoplus_hpca_srvy_l1b',...
  'Vhplus_dbcs_hpca_srvy_l1b','Vheplus_dbcs_hpca_srvy_l1b','Vheplusplus_dbcs_hpca_srvy_l1b','Voplus_dbcs_hpca_srvy_l1b',...
  'Phplus_dbcs_hpca_srvy_l1b','Pheplus_dbcs_hpca_srvy_l1b','Pheplusplus_dbcs_hpca_srvy_l1b','Poplus_dbcs_hpca_srvy_l1b',...
  'Thplus_dbcs_hpca_srvy_l1b','Theplus_dbcs_hpca_srvy_l1b','Theplusplus_dbcs_hpca_srvy_l1b','Toplus_dbcs_hpca_srvy_l1b',...
  'Vhplus_gsm_hpca_srvy_l1b','Vheplus_gsm_hpca_srvy_l1b','Vheplusplus_gsm_hpca_srvy_l1b','Voplus_gsm_hpca_srvy_l1b',...
  'Phplus_gsm_hpca_srvy_l1b','Pheplus_gsm_hpca_srvy_l1b','Pheplusplus_gsm_hpca_srvy_l1b','Poplus_gsm_hpca_srvy_l1b',...
  'Thplus_gsm_hpca_srvy_l1b','Theplus_gsm_hpca_srvy_l1b','Theplusplus_gsm_hpca_srvy_l1b','Toplus_gsm_hpca_srvy_l1b',...
  'Nhplus_hpca_brst_l2','Nheplus_hpca_brst_l2','Nheplusplus_hpca_brst_l2','Noplus_hpca_brst_l2',...
  'Tshplus_hpca_brst_l2','Tsheplus_hpca_brst_l2','Tsheplusplus_hpca_brst_l2','Tsoplus_hpca_brst_l2',...
  'Vhplus_dbcs_hpca_brst_l2','Vheplus_dbcs_hpca_brst_l2','Vheplusplus_dbcs_hpca_brst_l2','Voplus_dbcs_hpca_brst_l2',...
  'Phplus_dbcs_hpca_brst_l2','Pheplus_dbcs_hpca_brst_l2','Pheplusplus_dbcs_hpca_brst_l2','Poplus_dbcs_hpca_brst_l2',...
  'Thplus_dbcs_hpca_brst_l2','Theplus_dbcs_hpca_brst_l2','Theplusplus_dbcs_hpca_brst_l2','Toplus_dbcs_hpca_brst_l2',...
  'Vhplus_gsm_hpca_brst_l2','Vheplus_gsm_hpca_brst_l2','Vheplusplus_gsm_hpca_brst_l2','Voplus_gsm_hpca_brst_l2',...
  'Phplus_gsm_hpca_brst_l2','Pheplus_gsm_hpca_brst_l2','Pheplusplus_gsm_hpca_brst_l2','Poplus_gsm_hpca_brst_l2',...
  'Thplus_gsm_hpca_brst_l2','Theplus_gsm_hpca_brst_l2','Theplusplus_gsm_hpca_brst_l2','Toplus_gsm_hpca_brst_l2',...
  'Nhplus_hpca_brst_l1b','Nheplus_hpca_brst_l1b','Nheplusplus_hpca_brst_l1b','Noplus_hpca_brst_l1b',...
  'Tshplus_hpca_brst_l1b','Tsheplus_hpca_brst_l1b','Tsheplusplus_hpca_brst_l1b','Tsoplus_hpca_brst_l1b',...
  'Vhplus_dbcs_hpca_brst_l1b','Vheplus_dbcs_hpca_brst_l1b','Vheplusplus_dbcs_hpca_brst_l1b','Voplus_dbcs_hpca_brst_l1b',...
  'Phplus_dbcs_hpca_brst_l1b','Pheplus_dbcs_hpca_brst_l1b','Pheplusplus_dbcs_hpca_brst_l1b','Poplus_dbcs_hpca_brst_l1b',...
  'Thplus_dbcs_hpca_brst_l1b','Theplus_dbcs_hpca_brst_l1b','Theplusplus_dbcs_hpca_brst_l1b','Toplus_dbcs_hpca_brst_l1b',...
  'Vhplus_gsm_hpca_brst_l1b','Vheplus_gsm_hpca_brst_l1b','Vheplusplus_gsm_hpca_brst_l1b','Voplus_gsm_hpca_brst_l1b',...
  'Phplus_gsm_hpca_brst_l1b','Pheplus_gsm_hpca_brst_l1b','Pheplusplus_gsm_hpca_brst_l1b','Poplus_gsm_hpca_brst_l1b',...
  'Thplus_gsm_hpca_brst_l1b','Theplus_gsm_hpca_brst_l1b','Theplusplus_gsm_hpca_brst_l1b','Toplus_gsm_hpca_brst_l1b',...
  'Nhplus_hpca_srvy_sitl','Nheplus_hpca_srvy_sitl','Nheplusplus_hpca_srvy_sitl','Noplus_hpca_srvy_sitl',...
  'Tshplus_hpca_srvy_sitl','Tsheplus_hpca_srvy_sitl','Tsheplusplus_hpca_srvy_sitl','Tsoplus_hpca_srvy_sitl',...
  'Vhplus_dbcs_hpca_srvy_sitl','Vheplus_dbcs_hpca_srvy_sitl','Vheplusplus_dbcs_hpca_srvy_sitl','Voplus_dbcs_hpca_srvy_sitl',...
  'Phplus_dbcs_hpca_srvy_sitl','Pheplus_dbcs_hpca_srvy_sitl','Pheplusplus_dbcs_hpca_srvy_sitl','Poplus_dbcs_hpca_srvy_sitl',...
  'Thplus_dbcs_hpca_srvy_sitl','Theplus_dbcs_hpca_srvy_sitl','Theplusplus_dbcs_hpca_srvy_sitl','Toplus_dbcs_hpca_srvy_sitl',...
  'Vhplus_gsm_hpca_srvy_sitl','Vheplus_gsm_hpca_srvy_sitl','Vheplusplus_gsm_hpca_srvy_sitl','Voplus_gsm_hpca_srvy_sitl',...
  'Phplus_gsm_hpca_srvy_sitl','Pheplus_gsm_hpca_srvy_sitl','Pheplusplus_gsm_hpca_srvy_sitl','Poplus_gsm_hpca_srvy_sitl',...
  'Thplus_gsm_hpca_srvy_sitl','Theplus_gsm_hpca_srvy_sitl','Theplusplus_gsm_hpca_srvy_sitl','Toplus_gsm_hpca_srvy_sitl',...
  'Nhplus_hpca_sitl','aspoc_status',...
  'Omnifluxoplus_hpca_brst_l2','Omnifluxhplus_hpca_brst_l2','Omnifluxheplus_hpca_brst_l2','Omnifluxheplusplus_hpca_brst_l2',...
  ...'Omnifluxion_epd_feeps_brst_l2'...
  }; % XXX THESE MUST BE THE SAME VARS AS BELOW

if strcmp(varStr,'vars') % collect all vars, for testing
  res = vars;
  return
end

if isempty(intersect(varStr,vars))
  errS = ['variable not recognized: ' varStr];
  irf.log('critical',errS);
  vars = sort(vars);
  disp('Implemented vars are:')
  for iVar = 1:length(vars)
    fprintf('  %s\n',vars{iVar})
  end
  error(errS)
end

switch varStr
  case 'dfg_ql_srvy', varStr = 'B_dmpa_dfg_srvy_ql';
  case 'afg_ql_srvy', varStr ='B_dmpa_afg_srvy_ql';
  case 'Nhplus_hpca_sitl', varStr = 'Nhplus_hpca_srvy_sitl';
  case {'R_gse','R_gsm','V_gse','V_gsm'}
    vC = varStr(1); cS = varStr(3:5);
    
    if mmsId>0
      res = mms.db_get_ts(['mms' mmsIdS '_mec_srvy_l2_epht89d'],...
        ['mms' mmsIdS '_mec_' lower(vC) '_' cS],Tint);
      if ~isempty(res), return, end
      
      % LAST RESORT: Load position of MAG files
      if vC=='R'
      % Load from L2pre B
      res = mms.db_get_ts(...
        ['mms' mmsIdS '_dfg_srvy_l2pre'],['mms' mmsIdS '_pos_' cS],Tint);
      if ~isempty(res), return, end
      % Load from L2 B
      res = mms.db_get_ts(...
        ['mms' mmsIdS '_fgm_srvy_l2'],['mms' mmsIdS '_fgm_r_' cS '_srvy_l2'],Tint);
      if ~isempty(res), return, end
      % Load from QL B
      res = mms.db_get_ts(...
        ['mms' mmsIdS '_dfg_srvy_ql'],['mms' mmsIdS '_ql_pos_' cS],Tint);
      return
      end
    end
    
    % Do resampling similar to mms_update_ephemeris
    TintTmp = irf.tint(Tint.start+(-60),Tint.stop+60);
    TintTmp = EpochUnix([fix(TintTmp.start.epochUnix/60)*60 ...
      ceil(TintTmp.stop.epochUnix/60)*60]);
    TintTmp = EpochTT(TintTmp);
    res.time = EpochTT((TintTmp.start.epoch:int64(30*1e9):TintTmp.stop.epoch)');
    for mmsId=1:4
      mmsIdS = num2str(mmsId);
      
        dTmp = mms.db_get_ts(['mms' mmsIdS '_mec_srvy_l2_epht89d'],...
          ['mms' mmsIdS '_mec_' lower(vC) '_' cS],Tint);
      
      if isempty(dTmp) &&  vC=='V', continue, end
      
      if isempty(dTmp)
        % LAST RESORT: Load position of MAG files
        % Load from L2pre B
        dTmp = mms.db_get_ts(['mms' mmsIdS '_dfg_srvy_l2pre'],...
          ['mms' mmsIdS '_pos_' cS],TintTmp);
        if isempty(dTmp)
          % Load from L2 B
          dTmp = mms.db_get_ts(['mms' mmsIdS '_fgm_srvy_l2'],...
            ['mms' mmsIdS '_fgm_r_' cS  '_srvy_l2'],TintTmp);
        end
        if isempty(dTmp)
          % Load from QL B
          dTmp = mms.db_get_ts(['mms' mmsIdS '_dfg_srvy_ql'],...
            ['mms' mmsIdS '_ql_pos_' cS],TintTmp);
        end
      end
      if isempty(dTmp), continue, end
      dTmp = comb_ts(dTmp);
      dTmp.data = double(dTmp.data);
      dTmpR = dTmp.resample(res.time,'spline');
      res.([cS vC mmsIdS]) = dTmpR.data; 
    end 
    return
  case 'tetra_quality'
    % Begin looking for Def. quality
    quality = mms.db_get_variable('mms_ancillary_defq','quality',Tint);
    if isempty(quality)
      irf.log('warning', 'Did not find any definite tetrahedra quality. Looking for predicted.');
      list = mms.db_list_files('mms_ancillary_predq',Tint);
      if(~isempty(list))
        % Load the last predicted file to match Tint
        quality = mms_load_ancillary([list(end).path, filesep, list(end).name], 'predq');
      end
    end
    if(~isempty(quality))
      rTs = irf.ts_scalar(EpochTT(quality.time), quality.quality);
      res = rTs.tlim(Tint);
    end
    return
  case 'aspoc_status'
    dsetName = ['mms', mmsIdS, '_aspoc_srvy_l2'];
    pref = ['mms', mmsIdS, '_aspoc_status'];
    res = mms.db_get_ts(dsetName, pref, Tint);
    return
end

Vr = splitVs(varStr);
dsetName = ['mms' mmsIdS '_' Vr.inst '_' Vr.tmmode '_' Vr.lev];
compS = ''; pref = ''; suf = ''; pref_err = '';
err = [];

switch Vr.inst
  case 'fsm'
    vn = ['mms' mmsIdS '_' Vr.inst '_b_' Vr.cs '_' Vr.tmmode '_' Vr.lev];
    dsetName = [dsetName '_8khz'];
    res = mms.db_get_ts(dsetName, vn, Tint);
  case {'fgm','dfg','afg'}
    switch Vr.lev
      case 'l2'
        vn = ['mms' mmsIdS '_' Vr.inst '_b_' Vr.cs '_' Vr.tmmode '_' Vr.lev];
        res = mms.db_get_ts(dsetName, vn, Tint);
      case 'l2pre' 
        vn = ['mms' mmsIdS '_' Vr.inst '_b_' Vr.cs '_' Vr.tmmode '_' Vr.lev];
        res = mms.db_get_ts(dsetName, vn, Tint);
        % XXX: once there will be no v3x files, this entry can be combined
        % with 'l2'
        if isempty(res)
          irf.log('notice','Trying v3.x files')
          vn = ['mms' mmsIdS '_' Vr.inst '_' Vr.tmmode '_' Vr.lev '_' Vr.cs];
          res = mms.db_get_ts(dsetName, vn, Tint);
        end
      case 'ql'
        vn = ['mms' mmsIdS '_' Vr.inst '_' Vr.tmmode '_' Vr.cs];
        res = mms.db_get_ts(dsetName, vn, Tint);
      otherwise, error('should not be here')
    end
    if isempty(res), return, end
    if strcmp(Vr.lev,'srvy')
      ind = diff(res.time.ttns) <= 122000; % FIXME: what is brst min dt for A/DFG?
      if( sum(ind) < (length(res)-2) )
        % Remove samples that are too close, but ensure some output if only
        % two samples with very high sam  ple rate.
        irf.log('notice',['Removing ',sum(ind), ...
          ' samples due to overlap AFG/DFG when transitioning between fast/slow mode.']);
        res = res(~ind);
      end
    end
  case 'edi'
    switch Vr.param
      case 'E'
        pref = ['mms' mmsIdS '_edi_' lower(Vr.param) '_' Vr.cs '_' Vr.tmmode '_' Vr.lev];
        dset = 'efield';
        dsetName = ['mms' mmsIdS '_edi_' Vr.tmmode '_' Vr.lev '_' dset];
        res = mms.db_get_ts(dsetName,pref,Tint);
      case {'Flux-amb-pm2','Flux-err-amb-pm2'}
        err_pref = [];
        if strcmp(Vr.param,'Flux-err-amb-pm2')
          err_pref = '_delta';
          Vr.param = 'Flux-amb-pm2';
        end        
        % for amb-pm2, node 1 is closest to 0/180
        edi_tk = tokenize(Vr.param,'-');
        dsetName = [dsetName '_' edi_tk{2} '-' edi_tk{3}];
        nodes = 1:4;
        pitchangles = [0 180];
        flux = cell(numel(pitchangles),numel(nodes));
        for ipitchangle = 1:numel(pitchangles)
          for inode = 1:numel(nodes)
            pitchangle = pitchangles(ipitchangle);
            node = nodes(inode);
            pref = ['mms' mmsIdS '_' Vr.inst '_flux' num2str(inode)  '_' num2str(pitchangle) err_pref '_' Vr.tmmode '_' Vr.lev];
            flux{ipitchangle,inode} = get_ts('scalar');            
          end
        end
        if isempty(flux{1,1})
          res = [];          
          return
        end
        paddistarr = [flux{1,1}.data flux{1,2}.data flux{1,3}.data flux{1,4}.data flux{2,4}.data flux{2,3}.data flux{2,2}.data flux{2,1}.data];
        d_angle = 180/16;
        pitchangle_edges_edi = [0 1 2 3 4 12 13 14 15 16]*d_angle;
        pitchangle_centers_edi = [0.5 1.5 2.5 3.5 12.5 13.5 14.5 15.5]*d_angle;
        E_edi = 500;
        dE_edi = E_edi*0.1*0.5; % width is 10% of center energy 
        E_edges = E_edi + dE_edi*[-1 1];
        res = PDist(flux{1}.time,reshape(paddistarr,[size(paddistarr,1),1,size(paddistarr,2)]),'pitchangle',E_edi*ones(flux{1}.length,1),pitchangle_centers_edi);
        res.units = flux{1}.units;
        res.siConversion = flux{1}.siConversion;
        res.species = 'electrons';
        tmp_name = pref; tmp_name(14)='*'; tmp_name = strrep(tmp_name,'180','*'); tmp_name = strrep(tmp_name,'0','*');
        res.name = tmp_name;
        % ancillary data
        res.ancillary.dt_minus = 0.5*(flux{1}.time(2)-flux{1}.time(1));
        res.ancillary.dt_plus = 0.5*(flux{1}.time(2)-flux{1}.time(1));
        res.ancillary.energy = E_edi*ones(flux{1}.length,1);
        res.ancillary.energy0 = E_edi;
        res.ancillary.energy1 = E_edi;
        res.ancillary.esteptable = ones(flux{1}.length,1);
        res.ancillary.delta_energy_minus = dE_edi;
        res.ancillary.delta_energy_plus = dE_edi;
        res.ancillary.pitchangle_edges = pitchangle_edges_edi;  
        res.ancillary.delta_pitchangle_minus = d_angle*ones(1,8)*0.5;
        res.ancillary.delta_pitchangle_plus = d_angle*ones(1,8)*0.5;
      case 'Flux-amb-pm'
    end    
  case 'fpi'
    switch Vr.param(end)
      case 'i', sensor = 'dis';
      case 'e', sensor = 'des';
      otherwise 
        error('invalid specie')
    end
    
    switch Vr.lev
      case {'l2','l2pre','l1b'}
        if strcmp(Vr.param(1:2),'PD')
          dsetName = [dsetName '_' sensor '-dist'];
        else
            if length(Vr.param)>4
                if strcmp(Vr.param(1:4), 'part')
                    dsetName = [dsetName '_' sensor '-partmoms'];
                else
                    dsetName = [dsetName '_' sensor '-moms'];
                end
            else
                    dsetName = [dsetName '_' sensor '-moms'];
            end
        end
      case 'ql'
        dsetName = [dsetName '_' sensor];
      case 'sitl'
      otherwise, error('should not be here')
    end
    
    switch Vr.param
      case {'PDe','PDi', 'PDERRe', 'PDERRi'}
        switch Vr.lev
          case {'l2'}
            pref = ['mms' mmsIdS '_' sensor];        
          otherwise, error('should not be here')
        end        
        res = get_ts('skymap');
      case {'Ni','Ne'}
        switch Vr.lev
          case {'l2','l2pre'}
           % pref = ['mms' mmsIdS '_' sensor '_numberdensity_dbcs_' Vr.tmmode];
            % V3.1 FPI
            pref = ['mms' mmsIdS '_' sensor '_numberdensity_' Vr.tmmode];
            pref_err = ['mms' mmsIdS '_' sensor '_numberdensity_err_' Vr.tmmode];
          case 'l1b'
            pref = ['mms' mmsIdS '_' sensor '_numberdensity'];
          case 'ql'
            pref = ['mms' mmsIdS '_' sensor '_numberdensity_fast'];
          case 'sitl'
            pref = ['mms' mmsIdS '_fpi_' upper(sensor) 'numberDensity'];
          otherwise, error('should not be here')
        end
        res = get_ts('scalar');
      case {'partNi','partNe'}      % + on 20190909, wyli
        switch Vr.lev
          case {'l2'}
            % Only for l2 data
            pref = ['mms' mmsIdS '_' sensor '_numberdensity_part_' Vr.tmmode];
          otherwise, error('should not be here')
        end
        res = mms.db_get_ts(dsetName,pref,Tint);  % + on 20190909, wyli  
      case {'Enfluxi', 'Enfluxe'}
        switch Vr.lev
            case 'ql'
            pref = ['mms' mmsIdS '_' sensor '_energyspectr_omni_fast']; 
            otherwise, error('should not be here')
        end
        res = mms.db_get_ts(dsetName,pref,Tint);    
      case {'Energyi', 'Energye'}
        switch Vr.lev
            case 'ql'
            pref = ['mms' mmsIdS '_' sensor '_energy_fast']; 
            otherwise, error('should not be here')
        end
        res = mms.db_get_ts(dsetName,pref,Tint);            
      case {'Tsi','Tse'}
        getQ = 'trace';
        switch Vr.lev
          case {'l2','l2pre'}
            pref = ['mms' mmsIdS '_' sensor '_temp'];
            suf = ['_dbcs_' Vr.tmmode];
            compS = struct('xx','xx','yy','yy','zz','zz');
          case {'l1b','ql'}
            pref = ['mms' mmsIdS '_' sensor '_Temp'];
          case 'sitl'
            pref = ['mms' mmsIdS '_fpi_' upper(sensor) 'temp'];
            getQ = 'ts';
          otherwise, error('should not be here')
        end
        res = get_ts(getQ);
      case {'Tparai','Tparae', 'Tperpi','Tperpe'}
          tmpFAC = Vr.param(2:5);
        switch Vr.lev
          case {'l2'}
            pref = ['mms' mmsIdS '_' sensor '_temp' tmpFAC '_' Vr.tmmode];
          otherwise, error('should not be here')
        end
        res = get_ts('scalar');
      case {'partTparai','partTparae', 'partTperpi','partTperpe'} % + on 20190909
          tmpFAC = Vr.param(6:9);
        switch Vr.lev
          case {'l2'}
            pref = ['mms' mmsIdS '_' sensor '_temp' tmpFAC '_part_' Vr.tmmode];
            res = mms.db_get_ts(dsetName, pref,Tint);
            if ~isempty(res), return, end            
          otherwise, error('should not be here')
        end
      case {'Ti', 'Te', 'Pi', 'Pe'}
        % try to load v3.x
        switch Vr.param(1)
          case 'T', momType = 'temptensor'; % temperature
          case 'P', momType = 'prestensor'; % pressure
          otherwise, error('should not be here 2')
        end
        pref = ['mms' mmsIdS '_' sensor '_'];
        suf = ['_' Vr.cs '_' Vr.tmmode];
        res = mms.db_get_ts(dsetName,[pref momType suf],Tint);
        if ~isempty(res), return, end
        % continue with v2.x
        switch Vr.param(1)
          case 'T' % temperature
            momType = 'Temp';
          case 'P' % pressure
            momType = 'Pres';
          otherwise, error('should not be here 2')
        end
        switch Vr.lev
          case {'l2','l2pre'}
            pref = ['mms' mmsIdS '_' sensor '_' lower(momType)];
            suf = ['_' Vr.cs '_' Vr.tmmode];
            compS = struct('xx','xx','xy','xy','xz','xz','yy','yy','yz','yz','zz','zz'); 
          case {'l1b','ql'}
            pref = ['mms' mmsIdS '_' sensor '_' momType];
          case 'sitl'
            pref = ['mms' mmsIdS '_fpi_' upper(sensor) 'temp'];
          otherwise, error('should not be here')
        end
        res = get_ts('tensor2');
      case {'partTi', 'partTe', 'partPi', 'partPe'}     % + on 20190909
        % only load v3.x
        switch Vr.param(5)
          case 'T', momType = 'temptensor'; % temperature
          case 'P', momType = 'prestensor'; % pressure
          otherwise, error('should not be here 2')
        end
        pref = ['mms' mmsIdS '_' sensor '_'];
        suf = ['_part_' Vr.cs '_' Vr.tmmode];
        res = mms.db_get_ts(dsetName,[pref momType suf],Tint);
        if ~isempty(res), return, end  
        res = get_ts('tensor2');
      case {'Vi','Ve'}
        pref = ['mms' mmsIdS '_' sensor '_bulk'];
        switch Vr.lev
          case {'l2','l2pre'}
            suf = ['_' Vr.cs '_' Vr.tmmode];
            compS = struct('x','x','y','y','z','z');
            % try to load V3
            res = mms.db_get_ts(dsetName,[pref 'v' suf],Tint);
            if ~isempty(res), return, end
          case 'l1b'
          case 'ql'
            suf = ['_' Vr.cs '_' Vr.tmmode];
            compS = struct('x','x','y','y','z','z');
            % try to load V3
            res = mms.db_get_ts(dsetName,[pref 'v' suf],Tint);
            if ~isempty(res), return, end
          case 'sitl'
            pref = ['mms' mmsIdS '_fpi_' Vr.param(end) 'BulkV_'];
            suf = '_DSC';
          otherwise, error('should not be here')
        end
        res = get_ts('vector');
      case {'partVi','partVe'}                  % + on 20190909;
        pref = ['mms' mmsIdS '_' sensor '_bulk'];
        switch Vr.lev
          case {'l2'}
            suf = ['_part_' Vr.cs '_' Vr.tmmode];
            compS = struct('x','x','y','y','z','z');
            % try to load V3
            res = mms.db_get_ts(dsetName,[pref 'v' suf],Tint);
            if ~isempty(res), return, end
          otherwise, error('Only l2 partmoms avaiable now.')
        end        
      otherwise, error('should not be here')
    end    
  case 'hpca'
    % Some restructuring to include spectrograms
    % This could also be done in splitVs.
    all_ions = {'hplus','heplus','heplusplus','oplus'};
    ion_index = cellfun(@(s) ~isempty(strfind(Vr.param, s)), all_ions);
    ion = all_ions{find(ion_index)};    
    param = extractBefore(Vr.param,ion); 
    
    switch param
      case {'N','V','Ts','P','T'}
        dsetName = ['mms' mmsIdS '_hpca_' Vr.tmmode '_' Vr.lev '_moments'];
      case {'Omniflux'}
        dsetName = ['mms' mmsIdS '_hpca_' Vr.tmmode '_' Vr.lev '_ion'];
    end
    %param = Vr.param(1); ion = Vr.param(2:end); 
    %if ion(1)=='s', param = [param ion(1)]; ion = ion(2:end); end % Ts
    switch ion
      case {'hplus','heplus','heplusplus','oplus'}
      otherwise, error('unrecognized ion')
    end
    doPDist = 0; doOmni = 0;
    switch param
      case 'N', v = 'number_density';
      case 'V', v = 'ion_bulk_velocity';
      case 'Ts', v = 'scalar_temperature';
      case 'P', v = 'ion_pressure';
      case 'T', v = 'temperature_tensor';
      case 'Omniflux', v = 'flux'; doPDist = 1; doOmni = 1;
      otherwise, error('unrecognized param')
    end
    pref = ['mms' mmsIdS '_hpca_' ion '_' v];
    if Vr.to>0
      switch Vr.cs
        case 'gsm', pref = [pref '_GSM'];
        case 'dbcs'
        otherwise, error('invalid CS')
      end
    end
    
    if doPDist
      res = get_ts('skymap');      
    else
      res = mms.db_get_ts(dsetName,pref,Tint);
    end
    
    % XXX this needs to be investigated with HPCA
    if ~isempty(res)
      irf.log('warning','setting zeros to NaN for HPCA')
      res.data(res.data==0) = NaN;
      if (Vr.to>1)
        new_res = TSeries(res.time,res.data,'TensorOrder',Vr.to,'TensorBasis','xyz',...
           'repres',{'x','y','z'});
        new_res.coordinateSystem =  Vr.cs;    
        new_res.name = res.name;
        new_res.siConversion = res.siConversion;
        new_res.userData = res.userData;
        new_res.units = res.units;
        res = new_res;        
      end
    end
  case 'scm'
    switch Vr.lev
      case 'l2'
        otherwise
          error('not implemented yet')
    end
    dset = 'scb';
    param = 'acb';
    pref = ['mms' mmsIdS '_scm_' param '_' Vr.cs '_' dset '_' Vr.tmmode '_' Vr.lev];
    dsetName = ['mms' mmsIdS '_scm_' Vr.tmmode '_' Vr.lev '_' dset];
    res = mms.db_get_ts(dsetName, pref, Tint);
    res = comb_ts(res);
  case 'edp' 
    switch Vr.lev
      case 'sitl'
        switch Vr.param
          case 'E', dset = 'dce'; param = ['dce_xyz_' Vr.cs];
          case 'E2d', dset = 'dce2d'; param = ['dce_xyz_' Vr.cs];
          case 'V', dset = 'scpot';  param = 'scpot';
        end
        pref = ['mms' mmsIdS '_edp_' param '_' Vr.tmmode '_' Vr.lev];
      case 'ql'
        switch Vr.param
          case 'E', dset = 'dce'; param = ['dce_xyz_' Vr.cs];
          case 'E2d', dset = 'dce2d'; param = ['dce_xyz_' Vr.cs];
        end
        pref = ['mms' mmsIdS '_edp_' param];
      case 'l1b'
        switch Vr.param
          case 'E', dset = 'dce'; param = 'dce_sensor';
          case 'V', dset = 'dce'; param = 'dcv_sensor';
        end
        pref = ['mms' mmsIdS '_edp_' param];
      case 'l2a'
        dset = 'dce2d';
        switch Vr.param
          case 'Phase', param = 'phase';
          case {'Es12','Es34'}, param = ['espin_p' Vr.param(3:4)];
          case 'Adcoff', param = 'adc_offset';
          case {'Sdev12','Sdev34'}, param = ['sdevfit_p' Vr.param(5:6)];
          case 'E', param = 'dce'; 
        end
        pref = ['mms' mmsIdS '_edp_' param '_' Vr.tmmode '_' Vr.lev];
      otherwise
        switch Vr.param
          case 'E', dset = 'dce'; param = ['dce_' Vr.cs];
          case 'Epar', dset = 'dce'; param = 'dce_par_epar';
          case 'E2d', dset = 'dce2d'; param = ['dce_' Vr.cs]; 
          case 'V', dset = 'scpot'; param = 'scpot';
          case 'Vpsp', dset = 'scpot'; param = 'psp'; 
          case 'V6', dset = 'scpot';  param = 'dcv';
          otherwise, error('unrecognized param')
        end   
        pref = ['mms' mmsIdS '_edp_' param '_' Vr.tmmode '_' Vr.lev];
    end
    dsetName = ['mms' mmsIdS '_edp_' Vr.tmmode '_' Vr.lev '_' dset];
    res = mms.db_get_ts(dsetName,pref,Tint);
  case 'epd_feeps'
    all_species = {'ion','electron'};
    species_index = cellfun(@(s) ~isempty(strfind(Vr.param, s)), all_species);
    species = all_species{find(species_index)};    
    param = extractBefore(Vr.param,species); 
    
    % the files only has feeps, e.g.: mms1_feeps_brst_l2_ion_20170804093413_v6.1.2.cdf
    % but the variables have epd_feeps, e.g.: mms1_epd_feeps_brst_l2_ion_top_intensity_sensorid_6
    dsetName = ['mms' mmsIdS '_' extractAfter(Vr.inst,'_') '_' Vr.tmmode '_' Vr.lev '_' species];
    res = get_ts('feeps_omni');
  otherwise
    error('not implemented yet')
end

  function res = get_ts(dataType)
    res = [];
    switch dataType
      case 'scalar'
        rX = mms.db_get_ts(dsetName,pref,Tint);
        if isempty(rX)
          irf.log('warning',...
            ['No data for ' dsetName '(' pref ')'])
          return
        end
        rX = comb_ts(rX);
        res = irf.ts_scalar(rX.time, rX.data); % This removes all "userData"...
        res.name = [varStr '_' mmsIdS];
        res.units = rX.units;
        res.siConversion = rX.siConversion;
        rX_err = mms.db_get_ts(dsetName,pref_err,Tint);% error
        if isempty(rX_err)
          irf.log('warning',...
            ['No data for ' dsetName '(' pref_err ')'])
          return
        else
          rX_err = comb_ts(rX_err);
          err = rX_err.data;
        end
        res.userData.error = err;
      case 'vector'
        if isempty(compS), compS.x = 'X'; compS.y = 'Y'; compS.z = 'Z'; end
        rX = mms.db_get_ts(dsetName,[pref compS.x suf],Tint);
        if isempty(rX)
          irf.log('warning',...
            ['No data for ' dsetName '(' [pref compS.x suf] ')'])
          return
        end
        rX = comb_ts(rX);
        rY = comb_ts(mms.db_get_ts(dsetName,[pref compS.y suf],Tint));
        rZ = comb_ts(mms.db_get_ts(dsetName,[pref compS.z suf],Tint));
        res = irf.ts_vec_xyz(rX.time, [rX.data rY.data rZ.data]);
        res.coordinateSystem = Vr.cs;
        res.name = [varStr '_' mmsIdS];
        res.units = rX.units;
        res.siConversion = rX.siConversion;
      case 'trace'
        if isempty(compS), compS.xx = 'XX'; compS.yy = 'YY'; compS.zz = 'ZZ'; end
        rX = mms.db_get_ts(dsetName, [pref compS.xx suf],Tint); 
        if isempty(rX)
          irf.log('warning',...
            ['No data for ' dsetName '(' [pref compS.xx suf] ')'])
          return
        end
        rX = comb_ts(rX);
        rY = mms.db_get_ts(dsetName, [pref compS.yy suf],Tint); rY = comb_ts(rY);
        rZ = mms.db_get_ts(dsetName, [pref compS.zz suf],Tint); rZ = comb_ts(rZ);
        rX.data = (rX.data + rY.data + rZ.data)/3;
        res = irf.ts_scalar(rX.time, rX.data);
        res.name = [varStr '_' mmsIdS];
        res.units = rX.units;
        res.siConversion = rX.siConversion;
      case 'ts'
        if isempty(compS), compS.par = 'Para'; compS.perp = 'Perp'; end
        rX = mms.db_get_ts(dsetName, [pref compS.par suf],Tint); 
        if isempty(rX)
          irf.log('warning',...
            ['No data for ' dsetName '(' [pref compS.par suf] ')'])
          return
        end
        rX = comb_ts(rX);
        rY = mms.db_get_ts(dsetName, [pref compS.perp suf],Tint); rY = comb_ts(rY);
        rX.data = rX.data/3 + rY.data*2/3;
        res = irf.ts_scalar(rX.time, rX.data);
        res.name = [varStr '_' mmsIdS];
        res.units = rX.units;
        res.siConversion = rX.siConversion;
      case 'tensor2'
        if isempty(compS) 
          compS = struct('xx','XX','xy','XY','xz','XZ','yy','YY','yz','YZ','zz','ZZ'); 
        end
        rXX = mms.db_get_ts(dsetName,[pref compS.xx suf],Tint);
        if isempty(rXX),irf.log('warning',...
            ['No data for ' dsetName '(' [pref compS.xx suf] ')'])
          return
        end
        rXX = comb_ts(rXX);
        rXY = mms.db_get_ts(dsetName,[pref compS.xy suf],Tint);rXY = comb_ts(rXY);
        rXZ = mms.db_get_ts(dsetName,[pref compS.xz suf],Tint);rXZ = comb_ts(rXZ);
        rYY = mms.db_get_ts(dsetName,[pref compS.yy suf],Tint);rYY = comb_ts(rYY);
        rYZ = mms.db_get_ts(dsetName,[pref compS.yz suf],Tint);rYZ = comb_ts(rYZ);
        rZZ = mms.db_get_ts(dsetName,[pref compS.zz suf],Tint);rZZ = comb_ts(rZZ);
    
        rData = nan(rXX.length,3,3);
        rData(:,1,1) = rXX.data;
        rData(:,1,2) = rXY.data;
        rData(:,1,3) = rXZ.data;
        rData(:,2,1) = rXY.data;
        rData(:,2,2) = rYY.data;
        rData(:,2,3) = rYZ.data;
        rData(:,3,1) = rXZ.data;
        rData(:,3,2) = rYZ.data;
        rData(:,3,3) = rZZ.data;
        
        res = irf.ts_tensor_xyz(rXX.time, rData);
        res.name = [varStr '_' mmsIdS];
        res.units = rXX.units;
        res.siConversion = rXX.siConversion;
        res.coordinateSystem = Vr.cs;
      case 'skymap'
        switch [Vr.inst Vr.tmmode]
          case 'fpibrst'
            if (length(Vr.param) == 3)  
                dist = mms.db_get_variable(dsetName,[pref '_dist_' Vr.tmmode],Tint);
            elseif (length(Vr.param) == 6 && strcmp(Vr.param(3:5), 'ERR'))
                dist = mms.db_get_variable(dsetName,[pref '_disterr_' Vr.tmmode],Tint);
            end
            theta = dist.DEPEND_2.data;
            dist = mms.variable2ts(dist);
            dist = dist.tlim(Tint);
            energy_data = mms.db_get_variable(dsetName,[pref '_energy_' Vr.tmmode],Tint);   
            energy_delta_data = mms.db_get_variable(dsetName,[pref '_energy_delta_' Vr.tmmode],Tint); 
               
            % energy delta_minus/plus 
            % no 'DELTA_MINUS_VAR' or 'DELTA_PLUS_VAR' when loading
            % 'PDi_fpi_brst_l2' from mms.get_data; please let wenya know if
            % you change the following if ... else ... end section.
            % 2018-04-19.
            % CN: Problem with loading the delta +- from energy_data arises 
            % if Tint spans more than one burst file. mms.db_get_variable 
            % only loads one of them. See
            % size(energy_data.DELTA_MINUS_VAR.data)
            % size(energy_data.data)
            % Workaround is done by defaulating to variable 
            % 'energy_delta_data' for both plus and minus (energy is 
            % centered). If 'energy_delta_data' is empty, go back to 
            % previous DELTA_MINUS_VAR and DELTA_PLUS_VAR as before, still 
            % not good though. Problem should probably be addressed at a
            % lower level inside mms.db_get_variable.
            if not(isempty(energy_delta_data))
              denergy = mms.variable2ts(energy_delta_data); % delta_energy_plus == delta_energy_minus
              denergy = denergy.tlim(Tint); 
              energy_minus = denergy.data;
              energy_plus = denergy.data;
            elseif (isfield(energy_data, 'DELTA_MINUS_VAR') && isfield(energy_data, 'DELTA_PLUS_VAR'))              
              energy_minus = squeeze(energy_data.DELTA_MINUS_VAR.data);
              energy_plus = squeeze(energy_data.DELTA_PLUS_VAR.data);              
            else
                irf.log('warning','DELTA_MINUS_VAR/DELTA_PLUS_VAR is not loaded.')                
            end
            energy0 = mms.db_get_variable(dsetName,[pref '_energy0_' Vr.tmmode],Tint);
            energy1 = mms.db_get_variable(dsetName,[pref '_energy1_' Vr.tmmode],Tint);
            phi = mms.db_get_ts(dsetName,[pref '_phi_' Vr.tmmode],Tint);
            %theta = mms.db_get_variable(dsetName,[pref '_theta_' Vr.tmmode],Tint);
            stepTable = mms.db_get_ts(dsetName,[pref '_steptable_parity_' Vr.tmmode],Tint);
            if isempty(energy0)
              energymat = mms.db_get_ts(dsetName,[pref '_energy_' Vr.tmmode],Tint);
              if stepTable.data(1)
                energy1 = energymat.data(1,:);
                energy0 = energymat.data(2,:);
              else
                energy1 = energymat.data(2,:);
                energy0 = energymat.data(1,:);
              end
            else
              energy0 = energy0.data;
              energy1 = energy1.data;
            end
            res = irf.ts_skymap(dist.time,dist.data,[],phi.data,theta,'energy0',energy0,'energy1',energy1,'esteptable',stepTable.data);
            if (exist('energy_minus','var') && exist('energy_plus','var'))
                res.ancillary.delta_energy_minus = energy_minus;
                res.ancillary.delta_energy_plus = energy_plus;
            end
          case 'fpifast'
            %dist = mms.db_get_variable(dsetName,[pref '_dist_' Vr.tmmode],Tint);
            if (length(Vr.param) == 3)  
                dist = mms.db_get_variable(dsetName,[pref '_dist_' Vr.tmmode],Tint);
            elseif (length(Vr.param) == 6 && strcmp(Vr.param(3:5), 'ERR'))
                dist = mms.db_get_variable(dsetName,[pref '_disterr_' Vr.tmmode],Tint);
            end            
            phi = dist.DEPEND_1.data;
            theta = dist.DEPEND_2.data;
            dist = mms.variable2ts(dist);
            dist = dist.tlim(Tint);            
            energy = mms.db_get_ts(dsetName,[pref '_energy_' Vr.tmmode],Tint);            
            energy = energy.tlim(Tint);
            res = irf.ts_skymap(dist.time, dist.data, energy.data, phi, theta);
          case 'hpcabrst' % not implemented
            dist = mms.db_get_variable(dsetName,[pref],Tint);
            dist_ts = mms.variable2ts(dist);
            omni_data = squeeze(irf.nanmean(dist_ts.data,2));
            dist = PDist(dist_ts.time,omni_data,'omni',dist.DEPEND_2.data);
            res = dist;
            res.siConversion = dist_ts.siConversion;
            res.units = dist_ts.units;
            res.species = ion;
            % 
            % obj.depend{1} = args{1}; args(1) = []; obj.representation{1} = {'energy'};
        end
        switch [Vr.inst Vr.tmmode]
          case {'fpibrst','fpifast'}
            res.units = 's^3/cm^6';
            if strcmp(sensor(2),'e')          
              res.species = 'electrons';
            elseif strcmp(sensor(2),'i')
              res.species = 'ions';
            end            
            res.siConversion = '1e12';
          case {'hpcabrst'}
        end
          res.name = dsetName;
          res.userData = dist.userData;
      case 'feeps_omni'
        % FROM SPEDAS: Added by DLT on 31 Jan 2017: set unique energy and gain correction factors per spacecraft
        eEcorr = [14.0, -1.0, -3.0, -3.0]; % energy correction
        iEcorr = [0.0, 0.0, 0.0, 0.0]; % energy correction
        eGfact = [1.0, 1.0, 1.0, 1.0]; % gain factor
        iGfact = [0.84, 1.0, 1.0, 1.0]; % gain factor
        switch species
          case 'ion', sensors = 6:8;
            Ecorr = iEcorr;
            Gfact = iGfact;
          case 'electron'
            Ecorr = eEcorr;
            Gfact = eGfact;
            switch mode
              case 'brst', sensors = [1:5, 9:12];
              case 'fast', sensors = [3:5, 11:12];
              otherwise, error('invalid mode')
            end
          otherwise, error('invalid specie')
        end
        
        dsetPref= ['mms' mmsIdS '_' Vr.inst '_' Vr.tmmode '_' Vr.lev '_' species];
        
        eLow = mms.db_get_variable(dsetName, [species '_energy_lower_bound'],Tint);
        eUp = mms.db_get_variable(dsetName, [species '_energy_upper_bound'],Tint);
        energies = (eLow.data + eUp.data)/2. + eval([species(1) 'Ecorr(mmsId)']); % eV

        nSensors = length(sensors);
        for iSen = 1:nSensors
          sen = sensors(iSen); suf = sprintf('intensity_sensorid_%d',sen);
          %sen = sensors(iSen); suf = sprintf('count_rate_sensorid_%d',sen);
          sufMask = sprintf('sector_mask_sensorid_%d',sen);
          top = mms.db_get_ts(dsetName,[dsetPref '_top_' suf],Tint);
          mask = mms.db_get_ts(dsetName,[dsetPref '_top_' sufMask],Tint);
%           if all(size((mask.data))==size((top.data)))
%             top.data(logical(mask.data)) = NaN;
%           else
%             top.data(logical(repmat(mask.data,1,length(energies)))) = NaN; % obsolete?
%           end
          bot = mms.db_get_ts(dsetName,[dsetPref '_bottom_' suf],Tint);
          mask = mms.db_get_ts(dsetName,[dsetPref '_bottom_' sufMask],Tint);
%           if all(size((mask.data))==size((bot.data)))
%             bot.data(logical(mask.data)) = NaN;
%           else
%             bot.data(logical(repmat(mask.data,1,length(energies)))) = NaN;
%           end
          %c_eval([species(1) 'Tit?=top;'  species(1) 'Bit?=bot;'],sen) %
          %The usage of eval leads to error: Attempt to add "iTit6" to a static workspace.
          Tit{iSen} = top;
          Bit{iSen} = bot;
        end
        % omni
        %eval(['dTmp=' specie(1) 'Tit' num2str(sensors(1)) ';'])
        dataTmp = Tit{1};
        omniD = NaN( [size(dataTmp.data) nSensors*2]);
        for iSen = 1:nSensors
          omniD(:,:,iSen) = Tit{iSen}.data;
          omniD(:,:,nSensors+iSen) = Bit{iSen}.data;
          %c_eval(['omniD(:,:,iSen) = ' specie(1) 'Tit?.data;'...
          %  'omniD(:,:,nSensors+iSen) = ' specie(1) 'Bit?.data;'],sensors(iSen))
        end
        omnidata = dataTmp;
        omnidata.data = mean(double(omniD),3,'omitnan')*Gfact(mmsId);
        %eval([specie(1) 'Omni = dTmp; ' specie(1) 'Omni.data =' ...
        %  'mean(double(omniD),3,''omitnan'')*' specie(1) 'Gfact(ic);'])
                
        dist = PDist(omnidata.time,omnidata.data,'omni',energies*1e3); % energies keV -> eV
        res = dist;
        res.siConversion = Tit{1}.siConversion;
        res.units = Tit{1}.units;
        res.species = species;                    
      otherwise
        error('data type not implemented')
    end
  end
end %% MAIN

function TsOut = comb_ts(TsIn)
if ~iscell(TsIn), TsOut = TsIn; return; end
TsOut = TsIn{1};
for i=2:numel(TsIn)
  TsOut = combine(TsOut, TsIn{i});
end
end

function Res = splitVs(varStr)

tk = tokenize(varStr,'_');
nTk = length(tk);
if nTk <3 || nTk > 5, error('invalig STRING format'), end

hpcaParamsScal = {'Nhplus','Nheplus','Nheplusplus','Noplus',...
  'Tshplus','Tsheplus','Tsheplusplus','Tsoplus','Phase','Adcoff'}; 
hpcaParamSpec = {'Omnifluxoplus','Omnifluxhplus','Omnifluxheplus','Omnifluxheplusplus'};
hpcaParamsTens1 = {'Vhplus','Vheplus','Vheplusplus','Voplus'};
hpcaParamsTens2 = {'Phplus','Pheplus','Pheplusplus','Poplus',...
  'Thplus','Theplus','Theplusplus','Toplus'};
feepsParamsScal = {'Omnifluxion','Omnifluxelectron'};

  
param = tk{1};
switch param
  case {'Ni', 'partNi', 'Ne', 'partNe', 'Nhplus', 'Tsi', 'Tperpi', 'Tparai', 'partTperpi', 'partTparai', ...
          'Tse', 'Tperpe', 'Tparae', 'partTperpe', 'partTparae', 'PDe', 'PDi', 'PDERRe', 'PDERRi', 'V', 'V6', 'Vpsp', ...
      'Enfluxi', 'Enfluxe', 'Energyi', 'Energye', 'Epar', 'Sdev12', 'Sdev34','Flux-amb-pm2','Flux-err-amb-pm2'}
    tensorOrder = 0;
  case {'Vi', 'partVi', 'Ve', 'partVe', 'B', 'E','E2d','Es12','Es34'}
    tensorOrder = 1;
  case {'Pi', 'partPi', 'Pe', 'partPe', 'Ti', 'partTi', 'Te', 'partTe'}
    tensorOrder = 2;  
  case {hpcaParamsScal{:},hpcaParamSpec{:},feepsParamsScal{:}}
    tensorOrder = 0;
  case hpcaParamsTens1
    tensorOrder = 1;
  case hpcaParamsTens2
    tensorOrder = 2;
  otherwise 
    error('invalid PARAM')
end

coordinateSystem = []; idx = 1;
if tensorOrder > 0
  coordinateSystem = tk{idx+1}; idx = idx + 1;
  switch coordinateSystem
    case {'gse','gsm','dsl','dbcs','dmpa','ssc','bcs','par'}
    otherwise
      error('invalid COORDINATE_SYS')
  end
end

instrument = tk{idx+1}; idx = idx + 1;
switch instrument
  case {'fpi','edp','edi','hpca','fgm','dfg','afg','scm','fsm'}
  case {'epd'}
    instrument = [instrument '_' tk{idx+1}]; idx = idx + 1;    
  otherwise
    switch param
      case 'B', instrument = 'fgm'; idx = idx - 1;
      case {'Ve','Te','Ne','Pe'}, instrument = 'fpi'; idx = idx - 1;
      otherwise
        error('invalid INSTRUMENT')
    end
end

tmMode = tk{idx+1}; idx = idx + 1;
switch tmMode
  case {'brst', 'fast','slow','srvy'}
  otherwise
    tmMode = 'fast'; idx = idx - 1;
    irf.log('warning','assuming TM_MODE = FAST')
end

if length(tk)==idx, dataLevel = 'l2'; %default
else
  dataLevel = tk{idx+1};
  switch dataLevel
    case {'ql','sitl','l1b','l2a','l2pre','l2','l3'}
    otherwise
      error('invalid DATA_LEVEL level')
  end
end

Res = struct('param',param,'to',tensorOrder,'cs',coordinateSystem,...
  'inst',instrument,'tmmode',tmMode,'lev',dataLevel);
end