function [dataFixedPha,fixedPha,epochFixedPha] = mms_interp_fixed_pha(data, tEpoch2000, phaseDeg)

STEPS_PER_DEG = 1; phaShift=STEPS_PER_DEG/2;

epoch0 = tEpoch2000(1); epochTmp = double(tEpoch2000-epoch0);

if 1 % no gaps
  phaDegSeg = phaseDeg;
  epochSeg = epochTmp;
  dataSeg = data;
  
  [dataFixedPha,fixedPha,epochFixedPha] = interp_cont_segment();
  return
end

  function [dataFixedPhaSeg,phaFixedWrpSeg,epochFixedPhaSeg] = interp_cont_segment()
    phaRad = unwrap(phaDegSeg*pi/180);
    
    phaDegUnw = phaRad*180/pi;
    phaFixed = (fix(phaDegUnw(1)):STEPS_PER_DEG:fix(phaDegUnw(end)))' ...
      + phaShift;
    pha360 = 0:STEPS_PER_DEG:360; pha360 = pha360' + phaShift; pha360(end) = [];
    epochFixedPhaSeg = interp1(phaDegUnw,epochSeg,phaFixed,'linear');
    phaFixed(isnan(epochFixedPhaSeg)) = []; 
    epochFixedPhaSeg(isnan(epochFixedPhaSeg)) = [];
    phaFixedWrpSeg = mod(phaFixed,360);
    dataFixedPhaSeg = interp1(epochSeg,double(dataSeg),epochFixedPhaSeg,'spline');
  end
end