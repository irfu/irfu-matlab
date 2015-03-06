function quality = mms_sdp_bitmask2quality(dataType,bitmask)
% MMS_SDP_BITMASK2QUALITY  compute quality from bitmask
%
% QUALITY = MMS_SDP_BITMASK2QUALITY(DATA_TYPE, BITMASK)
%
% Compute QUALITY from BITMAS for DATA_TYPE ('e' or 'v')

switch lower(dataType)
  case 'e'
    %TODO: implement real functionality
    quality = bitmask*0;
  case 'v'
    %TODO: implement real functionality
    quality = bitmask*0;
  otherwise
    err_str = 'MMS_SDP_BITMASK2QUALITY invalid DATA_TYPE.';
    irf.log('critical', err_str);
    error('Matlab:MMS_SDP_BITMASK2QUALITY:INPUT', err_str);
end