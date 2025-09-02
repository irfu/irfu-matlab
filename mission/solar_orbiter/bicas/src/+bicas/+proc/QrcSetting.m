%
% Class that represents how to convert one particular QRC (QRCID) into
% modifications of quality ZVs, blanking of data etc. for an implicit group of
% DSIs. Implemented by subclasses which represent different groups of output
% DSIs. This implicitly means that the same QRC might apply to multiple
% instances (of different subclasses), e.g. for both L2 and L3 output.
%
% NOTE: The user of the subclasses should know their interface, and therefore
% this abstract superclass does not *NOT* define any interface (sic!). The
% reason is that different information is needed for different forms of
% processing.
% Ex: ZVs L2_QUALITY_BITMASK vs L3_QUALITY_BITMASK.
% Ex: ZV QUALITY_FLAG is technically shared between all output datasets, but
%     could be regarded as qualititatively different for different forms of
%     processing (i.e. warrants separate variables)
%
% NOTE: This superclass is still referenced by QRCSM by for assertions.
% NOTE: The class does not itself the QRCID or DSI themselves.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef(Abstract) QrcSetting
  % PROPOSAL: QRCS subclass naming which implies category of processing, not
  %           just output CDFs.
  %   NOTE: Cf. bicas.proc.*
  %     bicas.proc.L1L2/L2L2/L2L3
  %   L1L2
  %   L2L2 (for which there is currently no QRCS)
  %   L2L3
  %   L2L3Density



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(SetAccess=immutable)
    % *Cap* (max value) for the CDF ZV "QUALITY_FLAG" when the QRC applies.
    % NOTE: This is interpretation is in compliance with how the ZV
    % QUALITY_FLAG is supposed to be set/updated.
    qfl

    % GA "CAVEATS". Column array of strings.
    gaCaveats
  end



  methods(Access=public)



    function obj = QrcSetting(A)
      arguments
        A.qfl
        A.gaCaveats
      end

      assert(bicas.utils.validate_QFL(A.qfl))
      obj.qfl       =                 A.qfl;

      assert(iscolumn(A.gaCaveats) & isstring(A.gaCaveats))
      obj.gaCaveats = A.gaCaveats;
    end



  end



end
