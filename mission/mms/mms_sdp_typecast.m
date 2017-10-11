function [out] = mms_sdp_typecast(dataName, data)
% MMS_SDP_TYPECAST either get struct with classname as string in Matlab and
% in CDF representation along with fillvalue to be used for MMS "dataName".
% Or if "data" is supplied as second argument output is instead the "data"
% transformed to ensure it is of the correct class.
%
% Example of data conversion (transform what may be double (Matlab default)
% to proper class for DCE data (ie single)).
%    DCE_data = MMS_SDP_TYPECAST('dce', DATAC.dsl);
%
% Example to get the internal cdf classname (to be used when writing cdf
% files).
%    outStruct = MMS_SDP_TYPECAST('dce');
% outStruct has now the following fields:
%  .cdf     - class name in cdf representation (used when writing cdf)
%  .matlab  - class name in Matlab representation
%  .fillval - ISTP standard fillvalue for corresponding class (used when writing cdf)
%
% DataName: "epoch", "bitmask", "quality", "dce", "adc_offset", "scpot",
%  "phase", "delta_offset", "spinfits", "tensor_order", "deltaplus", 
%  "deltaminus" and "label".
%
% See also: MMS_SDP_CDFWRITE

% Matlab <=> CDF representation
% single <-> cdf_real4
% double <-> cdf_real8
% int8   <-> cdf_int1
% int16  <-> cdf_int2
% int32  <-> cdf_int4
% int64  <-> cdf_time_tt2000 (if it is EpochTT2000 times)
% int64  <-> cdf_int8        (if it is regular int64 data)
% uint8  <-> cdf_uint1
% uint16 <-> cdf_uint2
% uint32 <-> cdf_uint4
% char   <-> cdf_char (or cdf_uchar)

narginchk(1,2);
if(nargin==1), data=[]; end

switch(lower(dataName))
  case 'epoch'
    if(~isempty(data))
      out = int64(data);
    else
      out.cdf     = 'cdf_time_tt2000';
      out.matlab  = 'int64';
      out.fillval = int64(-9223372036854775808);
    end
  case {'deltaplus', 'deltaminus'}
    if(~isempty(data))
      out = int64(data);
    else
      out.cdf     = 'cdf_int8';
      out.matlab  = 'int64';
      out.fillval = int64(-9223372036854775808);
    end
  case 'bitmask'
    if(~isempty(data))
      out = uint16(data);
    else
      out.cdf     = 'cdf_uint2';
      out.matlab  = 'uint16';
      out.fillval = uint16(65535);
    end
  case 'quality'
    if(~isempty(data))
      out = uint8(data);
    else
      out.cdf     = 'cdf_uint1';
      out.matlab  = 'uint8';
      out.fillval = uint8(255);
    end
  case {'dce', 'adc_offset', 'scpot', 'phase', 'spinfits', 'delta_offset'}
    if(~isempty(data))
      out = single(data);
    else
      out.cdf     = 'cdf_real4';
      out.matlab  = 'single';
      out.fillval = single(-1.0E31);
    end
  case 'label'
    if(~isempty(data))
      errStr='Unexpected data input to label';
      irf.log('critical',errStr); error(errStr);
    end
    out.cdf    = 'cdf_char';
    out.matlab = 'char';
    % Labels dont use fillval
  case 'tensor_order'
    if(~isempty(data))
      out = int32(data);
    else
      out.cdf     = 'cdf_int4';
      out.matlab  = 'int32';
      out.fillval = int32(-2147483648);
    end
  otherwise
    % Should not be here
    errStr = ['Dataname ',dataName,' is not yet implemented.'];
    irf.log('critical', errStr);  error(errStr);
end

end
