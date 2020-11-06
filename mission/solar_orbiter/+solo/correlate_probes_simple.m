function OutTS = correlate_probes_simple(VDC)
%CORRELATE_PROBES_SIMPLE  Do linear fits betweent the probes signals


%%
Tstart = VDC.time.start+3600;
TstartS = Tstart.toUtc;
TstartS = [TstartS(1:11) '00:00:00Z'];
Tstart = EpochTT(TstartS);

Tstop = VDC.time.stop+3600;
TstopS = Tstop.toUtc;
TstopS = [TstopS(1:11) '00:00:00Z'];
Tstop = EpochTT(TstopS);

DT = 7200; % 2 hours
nSteps =  (Tstop-Tstart)/DT;

outTime = Tstart + ((1:nSteps)-0.5)*DT;

out = NaN(nSteps,7);
for i=1:nSteps
  Tint = irf.tint(Tstart + DT*(i-1),DT);
  Vtmp = VDC.tlim(Tint);
  v3 = double(Vtmp.z.data); % V3
  if isempty(v3) || length(v3)<1e4, continue, end
  idxNan = ~isnan(v3);
  v1 = double(Vtmp.x.data(idxNan)); v3 = v3(idxNan); % V1
  idxOk = abs(v3-median(v3))<5*std(v3);
  v2 = double(Vtmp.y.data(idxNan)); % V2
  d23 = median(v2) - median(v3); % V2 - V3
  v23 = (v3 + v2 -d23)/2; % V23
  out(i,1) = d23;
  d123 = median(v1) - median(v23); % V1 - V23
  out(i,2) = d123;
  v23corr = v23 + d123;
  % V23 vs V1
  k123 = v23corr(idxOk)\v1(idxOk);
  out(i,3) = k123;
end

OutTS = irf.ts_scalar(outTime,out);