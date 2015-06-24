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
  % Assume Good data to begin with.
  quality = 3*ones(size(bitmask));
  % Set quality to 1 (bad)
  bits_1 = MMS_CONST.Bitmask.BAD_BIAS + ...
    MMS_CONST.Bitmask.LOW_DENSITY_SATURATION + ...
    MMS_CONST.Bitmask.SWEEP_DATA + ...
    MMS_CONST.Bitmask.ADP_SHADOW;
  ind_1 = logical(bitand(bitmask, bits_1));
  quality(ind_1) = 1;
  % Set quality to 0 (really bad)
  bits_0 = bitor(MMS_CONST.Bitmask.SIGNAL_OFF, ...
    MMS_CONST.Bitmask.PROBE_SATURATION);
  ind_0 = logical(bitand(bitmask, bits_0));
  quality(ind_0) = 0;
end

end