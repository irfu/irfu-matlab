function dsl_off = mms_sdp_dsl_offset(scId, procId)
% Function to return offsets to be applied to DSL Ex (sunward) and DSL Ey
% (duskward) in processing of MMS despun electric fields.
% Example: 
%   MMS_CONST = mms_constants;
%   dsl_off = mms_sdp_dsl_offset(1, MMS_CONST.SDCProc.ql);
%  returns 
%   dsl_off.ex - Sunward offset (DSL Ex)
%   dsl_off.ey - Duskward offset (DSL Ey)
%
% Note: For now only a fixed value per s/c and for QL dce only.
%
% See also: MMS_SDP_DMGR

narginchk(2,2);
global MMS_CONST
if isempty(MMS_CONST), MMS_CONST = mms_constants(); end

switch procId
  case MMS_CONST.SDCProc.ql
    switch lower(scId)
      case 1
        dsl_off.ex = 2.2;
        dsl_off.ey = 0.0;
      case 2
        dsl_off.ex = 3.0;
        dsl_off.ey = 0.0;
      case 3
        dsl_off.ex = 2.55;
        dsl_off.ey = 0;
      case 4
        dsl_off.ex = 1.32;
        dsl_off.ey = 0;
      otherwise
        errStr = 'Unexpected scId';
        irf.log('critical',errStr); error(errStr);
    end
  otherwise
    irf.log('notice', 'ZERO sun-/duskward offsets in DSL E-field (sdp)');
    dsl_off.ex = 0.0;
    dsl_off.ey = 0.0;
end

end