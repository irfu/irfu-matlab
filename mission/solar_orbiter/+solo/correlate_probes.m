function OutTS =  correlate_probes(VDC,EDC,d23_inp,k23_inp, k123_inp, d123_inp)
%SOLO.CORRELATE_PROBES  Do LSQ fits betweent the probes signals
%
% OutTS = solo.correlate_probes(VDC,EDC,[d23],[k123])
%
% Inputs: VDC,EDC from L2 CWF files
%
% Optional inputs:
%   d23 - precomputed d23
%   k23 - precomputed k23
%  k123 - precomputed k123
%
% Output: Struct with the following TSeries
%    1 d23      : intercept of V2 vs V3 linear fit
%    2 k23      : slope of V2 vs V3 linear fit
%    3 d123     : median(V1) - median(V23);
%    4 k123     : k123 from LSQ, V1 = k123*V23 + del123
%    5 del123   : del123 from LSQ
%    6 d12      : median(V1-V2)=median(E12)


if nargin<3, flagD23 = false; else, flagD23 = true; end
if nargin<4, flagK23 = false; else, flagK23 = true; end
if nargin<5, flagK123 = false; else, flagK123 = true; end
if nargin<6, flagD123 = false; else, flagD123 = true; end
%%
Tstart = VDC.time.start+3600;
TstartS = Tstart.toUtc;
TstartS = [TstartS(1:11) '00:00:00Z'];
Tstart = EpochTT(TstartS);

Tstop = VDC.time.stop+3600;
TstopS = Tstop.toUtc;
TstopS = [TstopS(1:11) '00:00:00Z'];
Tstop = EpochTT(TstopS);

DT = 21600; % 6 hour sliding window. (+-3 hours)
nSteps =  (Tstop-Tstart)/DT;

outTime = Tstart + ((1:0.5:nSteps)-0.5)*DT;
out = NaN(outTime.length,6);

if flagD23
  d23R = d23_inp.resample(VDC);
  d23 = d23_inp.resample(outTime,'nearest'); out(:,1) = d23.data;
end
if flagK23
  k23R = k23_inp.resample(VDC);
  k23 = k23_inp.resample(outTime,'nearest'); out(:,6) = k23.data;
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
    %       d23 = median(e23);
    if ~isempty(v2) %if no data, set k23, d23 to NaN
      [p23,~]=polyfit(v2,v3,1);
      k23  = p23(1);
      d23 = p23(2);
    else
      k23=NaN;
      d23=NaN;
    end
    out(i,6) = k23;
    v23 = (v3 + k23*v2 +d23)/2; % V23 = (V2+V3)/2
    out(i,1) = d23;
  end

  V2_scaled = v2.*k23+d23; % Scale V2 to V3.
  V23 = (V2_scaled + v3)/2; % (V2 + V3) /2
  d123 = median(v1) - median(V23); % V1 - V23
  out(i,2) = d123;
  if flagK123
    K123 = k123R.tlim(Tint); k123 = K123.data(idxNan);
    del123 = median(v1) - median(V23.*k123); % V1 - V23

  else %Find the scaling between V23 and V1
    [k123,del123] = lsqfitgm(V23(idxOk),v1(idxOk));
    out(i,3) = k123;
  end
  out(i,4) = del123;

end

OutTS.d23 = irf.ts_scalar(outTime,out(:,1));
OutTS.d123 = irf.ts_scalar(outTime,out(:,2));
OutTS.k123 = irf.ts_scalar(outTime,out(:,3));
OutTS.del123 = irf.ts_scalar(outTime,out(:,4));
OutTS.d12 = irf.ts_scalar(outTime,out(:,5));
OutTS.k23 = irf.ts_scalar(outTime,out(:,6));

