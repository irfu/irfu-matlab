function OutTS = correlate_probes_simple(VDC,d23_inp, k123_inp)
%CORRELATE_PROBES_SIMPLE  Do linear fits betweent the probes signals


if nargin<2, flagD23 = false; else, flagD23 = true; end
if nargin<3, flagK123 = false; else, flagK123 = true; end

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
out = NaN(nSteps,4);

outTime = Tstart + ((1:nSteps)-0.5)*DT;
if flagD23
  d23R = d23_inp.resample(VDC);
  d23 = d23_inp.resample(outTime,'nearest'); out(:,1) = d23.data;
end
if flagK123
  k123R = k123_inp.resample(VDC);
  k123 = k123_inp.resample(outTime,'nearest'); out(:,3) = k123.data;
end


for i=1:nSteps
  Tint = irf.tint(Tstart + DT*(i-1),DT);
  Vtmp = VDC.tlim(Tint);
  v3 = double(Vtmp.z.data); % V3
  if isempty(v3) || length(v3)<1e4, continue, end
  idxNan = ~isnan(v3);
  v1 = double(Vtmp.x.data(idxNan)); v3 = v3(idxNan); % V1
  idxOk = abs(v3-median(v3))<5*std(v3);
  v2 = double(Vtmp.y.data(idxNan)); % V2
  if flagD23
    D23 = d23R.tlim(Tint); d23 = D23.data(idxNan);
    v23 = (v3 + v2 -d23)/2; % V23
  else
    d23 = median(v2) - median(v3); % V2 - V3
    v23 = (v3 + v2 -d23)/2; % V23
    out(i,1) = d23;
  end
  d123 = median(v1) - median(v23); % V1 - V23
  out(i,2) = d123;
  if flagK123
    K123 = k123R.tlim(Tint); k123 = K123.data(idxNan);
    del123 = median(v1) - median(v23.*k123); % V1 - V23
  else
    [k123,del123] = lsqfitgm(v23(idxOk),v1(idxOk));
    %v23corr = v23 + d123;
    % V23 vs V1
    %k123 = v23corr(idxOk)\v1(idxOk);
    out(i,3) = k123;
  end
  out(i,4) = del123;
end

OutTS = irf.ts_scalar(outTime,out);