%
% QRCS for producing L3 datasets (not just density).
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef QrcSettingL3 < bicas.proc.QrcSetting
  % PROPOSAL: Add QUALITY_FLAG?
  % NOTE: L3_QUALITY_BITMASK is only in L3 density.



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(SetAccess=immutable)
    % Column arrays of indices in the second dimension. Specifies which
    % components should be blanked.

    % For blanking L2 ZV VDC/EDC before being used before data is passed on to
    % solo.vdccal() and psp2ne().
    % --
    % NOTE: L2-->L3 processing does not use AC data.
    vdcFvIndexAr
    edcFvIndexAr

    % For blanking EFIELD data (after derivation; not intput).
    efieldFvIndexAr
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    function obj = QrcSettingL3(A)
      arguments
        A.qfl             = bicas.const.qrc.QFL_MAX
        A.vdcFvIndexAr    = zeros(0, 1);
        A.edcFvIndexAr    = zeros(0, 1);
        A.efieldFvIndexAr = zeros(0, 1);
        A.gaCaveats       = string.empty(0, 1);
      end

      obj@bicas.proc.QrcSetting(qfl=A.qfl, gaCaveats=A.gaCaveats);

      obj.assert_fvIndexAr(A.vdcFvIndexAr, 3)
      obj.vdcFvIndexAr    = A.vdcFvIndexAr;

      obj.assert_fvIndexAr(A.edcFvIndexAr, 3)
      obj.edcFvIndexAr    = A.edcFvIndexAr;

      obj.assert_fvIndexAr(A.efieldFvIndexAr, 3)
      obj.efieldFvIndexAr = A.efieldFvIndexAr;
    end



  end



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    function assert_fvIndexAr(fvIndexAr, nIndices)
      assert(iscolumn(      fvIndexAr))
      assert(all(ismember(  fvIndexAr, 1:nIndices)))
      irf.assert.number_set(fvIndexAr)
    end



  end    % methods(Static)



end
