function ADC_off = mms_sdp_adc_off(time, spinfits, ~)
%MMS_SDP_ADC_OFF  Compute ADC (raw data) offsets
%
% Compute ADC offset for each time stamp in DCE from spinfits
% See also: ClusterProc.getData, irf_waverage

% default settings
%nPointsADCOffset = 5; %or 7 or 9 or...?;
nPointsADCOffset = 21; % MMS

narginchk(3,3);
if isempty(time) || mms_is_error(time)
  errStr='Bad TIME input, cannot proceed.';
  irf.log('critical',errStr); error(errStr);
end
if isempty(spinfits) || mms_is_error(spinfits)
  errStr='Bad SPINFITS input, cannot proceed.';
  irf.log('critical',errStr); error(errStr);
end
if isempty(spinfits.time)
  % Empty spinfits could be caused by to short time series (Burst less than
  % 5 seconds processed without L2a dce2d ready, ie "QL"). One example
  % is mms2_edp_brst_l1b_dce_20161205125604_v1.4.0.cdf (size:344K, 6080
  % reconds ie less than 1 second duration).
  errStr='Empty spinfits, cannot proceed.';
  irf.log('critical', errStr); error(errStr);
end

sdpProbes = fieldnames(spinfits.sfit); % default {'e12', 'e34'}

for iProbe=1:numel(sdpProbes)
  % adc_off = ["sfit timestamp", "sfit A-coeff"], where timestamp are by
  % default every 5 seconds (tt2000 int64). Convert both to "double" for
  % interp1 and similar things to work.
  adc_off = [double(spinfits.time), double(spinfits.sfit.(sdpProbes{iProbe})(:,1))];

  adc_off(:,2) = movmedian(adc_off(:,2), nPointsADCOffset,'omitnan');

  if(size(adc_off,1)==1)
    % Only one good adc_offset (possibly because only one spinfit).
    ADC_off.(sdpProbes{iProbe})(1:length(time),1) = adc_off(:,2);
  else
    % Resample adc offset to match up with dce timestamps
    ADC_off.(sdpProbes{iProbe}) = interp1(adc_off(:,1), adc_off(:,2), ...
      double(time), 'linear', 'extrap');
  end
end

end