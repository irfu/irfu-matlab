function [dataFixedPha,fixedPha,epoch2000FixedPha] = mms_interp_fixed_pha(data, tEpoch2000, phaseDeg)

STEPS_PER_DEG = 0.5; phaShift = 0; %phaShift=STEPS_PER_DEG/2;
PHASE_OFF = 30; % offset the phase so that the boom is 45 deg to the sun

epoch0 = tEpoch2000(1); epochTmp = double(tEpoch2000-epoch0);

if max(diff(epochTmp))>1e9*60/3.1/2 % Gap > 1/ spin period
  irf.log('critical','long data gap')
  error('long data gap')
else % no gaps
  phaDegSeg = phaseDeg;
  epochSeg = epochTmp;
  dataSeg = data;
  
  [dataFixedPha,fixedPha,epochFixedPha] = interp_cont_segment();
  epoch2000FixedPha = epoch0 + int64(epochFixedPha);
  return
end

  function [dataFixedPhaSeg,phaFixedWrpSeg,epochFixedPhaSeg] = interp_cont_segment()
    phaDegUnw = unwrap(phaDegSeg*pi/180)*180/pi + PHASE_OFF;
    %phaFixed = (fix(phaDegUnw(1)):STEPS_PER_DEG:fix(phaDegUnw(end)))' ...
    %  + phaShift;
    phaFixed = ((phaDegUnw(1)-rem(phaDegUnw(1),360)):STEPS_PER_DEG:...
      (phaDegUnw(end)-rem(phaDegUnw(end),360)+360))' + phaShift;
    epochFixedPhaSeg = interp1(phaDegUnw,epochSeg,phaFixed,'linear');
    %phaFixed(isnan(epochFixedPhaSeg)) = []; 
    %epochFixedPhaSeg(isnan(epochFixedPhaSeg)) = [];
    phaFixedWrpSeg = mod(phaFixed,360);
    dataFixedPhaSeg = interp1(epochSeg,double(dataSeg),epochFixedPhaSeg,'spline');
  end
end