function quality = mms_sdp_bitmask2quality(dataType,bitmask)
% MMS_SDP_BITMASK2QUALITY  compute quality from bitmask
%
% QUALITY = MMS_SDP_BITMASK2QUALITY(DATA_TYPE, BITMASK)
%
% Compute QUALITY from BITMAS for DATA_TYPE ('e' or 'v')

global MMS_CONST
if isempty(MMS_CONST), MMS_CONST = mms_constants(); end

switch lower(dataType)
  case 'e'
    %TODO: implement real functionality
    quality = comp_quality(bitmask);
  case 'v'
    %TODO: implement real functionality
    quality = comp_quality(bitmask);
  otherwise
    err_str = 'MMS_SDP_BITMASK2QUALITY invalid DATA_TYPE.';
    irf.log('critical', err_str);
    error('Matlab:MMS_SDP_BITMASK2QUALITY:INPUT', err_str);
end

function quality = comp_quality(bitmask)
  % Assume Good data to begin with. And use the same type as bitmask at
  % first (for bitand comparison to work).
  quality = 3*ones(size(bitmask), ...
    getfield(mms_sdp_typecast('bitmask'),'matlab'));
  % Set quality to 2 (caution)
  bits_2 = MMS_CONST.Bitmask.ASPOC_RUNNING + ...
    MMS_CONST.Bitmask.ASYMM_CONF + ...
    MMS_CONST.Bitmask.BAD_BIAS + ...
    MMS_CONST.Bitmask.SW_WAKE_REMOVED;
  ind_2 = logical(bitand(bitmask, bits_2));
  quality(ind_2) = 2;
  % Set quality to 1 (bad)
  bits_1 = MMS_CONST.Bitmask.LOW_DENSITY_SATURATION + ...
    MMS_CONST.Bitmask.SWEEP_DATA + ...
    MMS_CONST.Bitmask.ADP_SHADOW;
  ind_1 = logical(bitand(bitmask, bits_1));
  quality(ind_1) = 1;
  % Set quality to 0 (really bad)
  bits_0 = MMS_CONST.Bitmask.SIGNAL_OFF + ...
    MMS_CONST.Bitmask.PROBE_SATURATION + ...
    MMS_CONST.Bitmask.ECLIPSE + ...
    MMS_CONST.Bitmask.MANEUVERS;
  ind_0 = logical(bitand(bitmask, bits_0));
  quality(ind_0) = 0;
  % Recast to proper quality datatype.
  quality = mms_sdp_typecast('quality', quality);
end

end
