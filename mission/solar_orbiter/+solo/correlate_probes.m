function OutTS =  correlate_probes(VDC,EDC,d23_inp,k23_inp,Gamma_inp, k123_inp, d123_inp)
%SOLO.CORRELATE_PROBES  Do LSQ fits betweent the probes signals
%
% OutTS = solo.correlate_probes(VDC,EDC,[d23],[k123])
%
% Inputs: VDC,EDC from L2 CWF files
%
% Optional inputs: 
%   d23 - precomputed d23
%  k123 - precomputed k123
%
% Output TSeries with the following columns:
%    1 d23    : median(V2-V3)=median(E23)
%    2 d123   : median(V1) - median(V23);
%    3 k123   : k123 from LSQ, V1 = k123*V23 + del123
%    4 del123 : del123 from LSQ
%    5 d12    : median(V1-V2)=median(E12)


if nargin<3, flagD23 = false; else, flagD23 = true; end
if nargin<4, flagK23 = false; else, flagK23 = true; end
if nargin<5, flagGamma = false; else, flagGamma = true; end
if nargin<6, flagK123 = false; else, flagK123 = true; end
if nargin<7, flagD123 = false; else, flagD123 = true; end
%%
Tstart = VDC.time.start+3600;
TstartS = Tstart.toUtc;
TstartS = [TstartS(1:11) '00:00:00Z'];
Tstart = EpochTT(TstartS);

Tstop = VDC.time.stop+3600;
TstopS = Tstop.toUtc;
TstopS = [TstopS(1:11) '00:00:00Z'];
Tstop = EpochTT(TstopS);

DT = 7200*3; % 2 hours
nSteps =  (Tstop-Tstart)/DT;

outTime = Tstart + ((1:0.5:nSteps)-0.5)*DT;
out = NaN(outTime.length,9);

cc_crit=0.7; %Only remove common mode when the correlation is greater than cc_crit.

if flagD23
  d23R = d23_inp.resample(VDC);
  d23 = d23_inp.resample(outTime,'nearest'); out(:,1) = d23.data;
end
if flagK23
  k23R = k23_inp.resample(VDC);
  k23 = k23_inp.resample(outTime,'nearest'); out(:,9) = k23.data;    
end
if flagGamma
  Gamma0R = Gamma_inp.Gamma0.resample(VDC);
  Gamma0 = Gamma_inp.Gamma0.resample(outTime,'nearest'); out(:,7)=Gamma0.data;
  Gamma1R = Gamma_inp.Gamma1.resample(VDC);
  Gamma1 = Gamma_inp.Gamma1.resample(outTime,'nearest'); out(:,6)=Gamma1.data;
  CCR = Gamma_inp.cc.resample(VDC);
  CC = Gamma_inp.cc.resample(outTime,'nearest'); out(:,8)=CC.data;
end
if flagK123
  k123R = k123_inp.resample(VDC);
  k123 = k123_inp.resample(outTime,'nearest'); out(:,3) = k123.data;
end
if flagD123
    d123R=d123_inp.resample(VDC);
    d123 = d123_inp.resample(outTime,'nearest');
end


for i=1:outTime.length
  Tint = irf.tint(outTime(i)+(-DT/2),DT);
  Vtmp = VDC.tlim(Tint); Etmp = EDC.tlim(Tint);
  v3 = double(Vtmp.z.data); % V3
  if isempty(v3) || length(v3)<1e4, continue, end
  idxNan = ~isnan(v3);
  v1 = double(Vtmp.x.data(idxNan)); v3 = v3(idxNan); % V1
  idxOk = abs(v3-median(v3))<5*std(v3);
  v2 = double(Vtmp.y.data(idxNan)); % V2
  e12 = double(Etmp.x.data(idxNan)); % E12
  e23 = double(Etmp.z.data(idxNan)); % E23
  
  d12 = median(e12); out(i,5) = d12; 
  if flagD23
    D23 = d23R.tlim(Tint); d23 = D23.data(idxNan);
    K23 = k23R.tlim(Tint); k23 = K23.data(idxNan);
    v23 = (v3 + v2 -d23)/2; % V23
  else
      d23 = median(e23);
      [p23,~]=polyfit(v2,v3,1);
      k23  = p23(1);
      out(i,9) = k23;
      v23 = (v3 + v2 -d23)/2; % V23
      out(i,1) = d23;
  end
  if flagGamma
      GAMMA1=Gamma1R.tlim(Tint); Gamma1=GAMMA1.data(idxNan);
      GAMMA0=Gamma0R.tlim(Tint); Gamma0=GAMMA0.data(idxNan);
      CC=CCR.tlim(Tint); cc=CC.data(idxNan);
      V2corr = v2-double(d23);
      V23 = (V2corr+v3)/2;
      V2cmr = V2corr-(Gamma0+Gamma1.*V23)/2;
      V3cmr = v3+(Gamma0+Gamma1.*V23)/2;
      v23cmr = (V3cmr+V2cmr)/2;
  else
      V2corr = v2 -double(d23);
      V23 = (V2corr + v3)/2; % (V2 + V3) /2
      V_delta23 = V2corr-v3;
      [p,s]=polyfit(V23,V_delta23,1); %p(1)=slope. p(2)=offset.
      cc=corrcoef(V23,V_delta23);
      % Now we know gamma. 
      if abs(cc)>cc_crit
          V2cmr = V2corr-(p(1)*V23+p(2))/2; %cmr = CommonModeRemoved
          V3cmr = v3+(p(1)*V23+p(2))/2; %cmr = CommonModeRemoved
          
          v23cmr = (V3cmr+V2cmr)/2; %V23 new.
      else
          p(1)=0;
          p(2)=0;
          v23cmr = v23; %V23 new.
      end
      
      out(i,6)=p(1); % V23 E23 slope
      out(i,7)=p(2); % V23 E23 offset
      out(i,8)=cc(1,2); %V23 E23 correlation coefficient.
  end

  d123 = median(v1) - median(v23cmr); % V1 - V23
  out(i,2) = d123;
  if flagK123
    K123 = k123R.tlim(Tint); k123 = K123.data(idxNan);
    del123 = median(v1) - median(v23cmr.*k123); % V1 - V23
    
  else
    [k123,del123] = lsqfitgm(v23cmr(idxOk),v1(idxOk));
    %v23corr = v23 + d123;
    % V23 vs V1
    %k123 = v23corr(idxOk)\v1(idxOk);
    out(i,3) = k123;
  end
  out(i,4) = del123;

end

OutTS.d23 = irf.ts_scalar(outTime,out(:,1));
OutTS.k23 = irf.ts_scalar(outTime,out(:,9));
OutTS.d123 = irf.ts_scalar(outTime,out(:,2));
OutTS.k123 = irf.ts_scalar(outTime,out(:,3));
OutTS.del123 = irf.ts_scalar(outTime,out(:,4));
OutTS.d12 = irf.ts_scalar(outTime,out(:,5));
OutTS.PotSlope = irf.ts_scalar(outTime,out(:,6));
OutTS.PotOffset = irf.ts_scalar(outTime,out(:,7));
OutTS.PotCorrcoeff = irf.ts_scalar(outTime,out(:,8));


