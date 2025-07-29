%
% QRCS for producing L3 datasets (not just density).
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef QrcSettingL3 < bicas.proc.QrcSetting
  % PROPOSAL: Add QUALITY_FLAG?
  % NOTE: L3_QUALITY_BITMASK is only in L3 density.



  properties(SetAccess=immutable)
    % Column arrays of column indices into the L2 input ZV VDC/EDC (i.e. in
    % interval 1-3) which are used for deriving L3 data. Specifies which
    % components should be blanked before the information is passed on to
    % solo.vdccal() and psp2ne().
    vdcFvIndexAr
    edcFvIndexAr
  end



  methods(Access=public)



    function obj = QrcSettingL3(A)
      arguments
        A.vdcFvIndexAr = zeros(0, 1);
        A.edcFvIndexAr = zeros(0, 1);
      end

      assert(iscolumn(      A.vdcFvIndexAr))
      assert(all(ismember(  A.vdcFvIndexAr, [1, 2, 3])))
      irf.assert.number_set(A.vdcFvIndexAr)
      obj.vdcFvIndexAr =    A.vdcFvIndexAr;

      assert(iscolumn(      A.edcFvIndexAr))
      assert(all(ismember(  A.edcFvIndexAr, [1, 2, 3])))
      irf.assert.number_set(A.edcFvIndexAr)
      obj.edcFvIndexAr =    A.edcFvIndexAr;
    end



  end



end
