function ADC_off = mms_sdp_adc_off(time,spinfits)
%MMS_SDP_ADC_OFF  Compute ADC (raw data) offsets
%
% Compute ADC offset for each time stamp in DCE from spinfits
% See also: ClusterProc.getData, irf_waverage

% default settings
FLAG_ADC_OFF_DESPIKE = true;
%nPointsADCOffset = 5; %or 7 or 9 or...?;
nPointsADCOffset = 21; % MMS

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

  max_off = adc_off(~isnan(adc_off(:,2)),:);
  adc_off_mean = mean(max_off(:,2));
  
  % Replace NaN with mean value
  adc_off(isnan(adc_off(:,2)),2) = adc_off_mean;

  if(FLAG_ADC_OFF_DESPIKE)
    max_off = 3*std(max_off(:,2));
    % if adc_despike, locate large adc_off
    idx = find( abs( adc_off(:,2) - adc_off_mean ) > max_off );
    if(~isempty(idx))
      adc_off(idx, 2) = 0;
      adc_off_mean = mean( adc_off( abs( adc_off(:, 2) )>0, 2) );
      adc_off(idx, 2) = adc_off_mean;
    end
  end

  % Smooth ADC offsets
  adc_off = mms_wavege(adc_off, nPointsADCOffset);
  
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


function [ out ] = mms_wavege(data, nPoints)
  % Weigted average function.
  narginchk(2,2); nargoutchk(1,1);

  if( ~ismember(nPoints,[5 7 9 21]) )
    errStr='nPoints must be 5, 7, 9 or 21';
    irf.log('critical',errStr);  error(errStr);
  end
  if( size(data,1)<=1 )
    irf.log('warning',['Not enough points (' num2str(size(data,1)) ') to average.']);
    out = data;
    return
  end

  ndata = size(data,1);
  ncol = size(data,2);

  out = data;
  out(isnan(out)) = 0; % set NaNs to zeros

  padd = zeros(1, floor(nPoints/2)); % Calculate padding
  for col=2:ncol
    dtmp = [padd, out(:,col)', padd]; % Apply padding at begining and end
    for j=1:ndata
      out(j,col) = w_ave(dtmp(j:j+nPoints-1), nPoints);
    end
  end

  % Make sure we do return matrix of the same size
  out = out(1:ndata, :);
end

function av = w_ave(x, nPoints)
  switch nPoints
    % get weight factor m, (normalized to one).
    case 5
      % Weights based on Cluster EFW
      m = [.1 .25 .3 .25 .1];
    case 7
      % Weights bases on Cluster EFW
      m = [.07 .15 .18 .2 .18 .15 .07];
    case 9
      % Weights based on almost "Binominal" distribution, given by:
      %y=binopdf(0:10,10,0.5);y(5)=y(5)+y(1);y(1)=0;y(7)=y(7)+y(end);y(end)=0;
      %m=y(2:end-1); sum(m)==1;
      m = [0.009765625, 0.0439453125, 0.1181640625, 0.205078125, ...
        0.2470703125, 0.205078125 0.1171875 0.0439453125 0.009765625];
    case 21
      % MMS "21 spinfits, made every 5 seconds" correspond to around 5
      % compelete spin revolutions.
      % Weights based on almost normal distribution, with rounded and
      % with slighlty more weight given to beginning and end.
      % Ensuring sum(m)==1
      m = [...
        0.012, 0.02, 0.024, 0.03, 0.04, 0.05, 0.06, 0.07, 0.075, 0.079, ...
        0.08, ...
        0.079, 0.075, 0.07, 0.06, 0.05, 0.04, 0.03, 0.024, 0.02, 0.012];
    otherwise
      errStr='nPoints must be 5, 7, 9 or 21';
      irf.log('critical',errStr);      error(errStr);
  end

  cor = sum(m(x==0)); % find missing points==0
  if cor==1
    av = 0;
  else
    av = sum(x.*m)/(1-cor);
  end

end
