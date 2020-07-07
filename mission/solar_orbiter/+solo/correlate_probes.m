function OutTS = correlate_probes(fName,d23)
%CORRELATE_PROBES  Do linear fits betweent the probes signals

if nargin<2, d23 = []; end

d = dataobj(fName);
VDC = get_ts(d,'VDC');

%%
Tstart = VDC.time(1)+3600;
TstartS = Tstart.toUtc;
TstartS = [TstartS(1:11) '00:00:00Z'];
Tstart = EpochTT(TstartS);

outTime = Tstart + ((1:12)-0.5)*7200;
if isempty(d23) % all 3 probes
  out = NaN(12,7);
  for i=1:12
    Tint = irf.tint(Tstart + 7200*(i-1),7200);
    Vtmp = VDC.tlim(Tint);
    x = double(Vtmp.z.data);
    if isempty(x) || length(x)<1e4, continue, end
    idxNan = ~isnan(x);
    y = double(Vtmp.x.data(idxNan)); x = x(idxNan);
    idxOk = abs(x-median(x))<5*std(x);
    p = polyfit(x(idxOk),y(idxOk),1);
    out(i,1:2) = p;
    %p = polyfit(x(~isnan(x)),Vtmp.y.data(~isnan(x)),1);
    y2 = double(Vtmp.y.data(idxNan));
    p = polyfit(x(idxOk),y2(idxOk),1);
    out(i,3:4) = p;
    d23 = median(y2) - median(x);
    x = (x + y2 -d23)/2;
    out(i,5) = d23;
    p = polyfit(x(idxOk),y(idxOk),1);
    out(i,6:7) = p;
  end
else % correlate V1 to V2+V3
  out = NaN(12,2);
  for i=1:12
    Tint = irf.tint(Tstart + 7200*(i-1),7200);
    Vtmp = VDC.tlim(Tint);
    x = double(Vtmp.z.data);
    if isempty(x) || length(x)<1e4, continue, end
    idxNan = ~isnan(x);
    y = double(Vtmp.x.data(idxNan)); x = x(idxNan); 
    y2 = double(Vtmp.y.data(idxNan));

    d23R = d23.resample(Vtmp);
    x = (x + y2 -d23R.data(idxNan))/2; % (V2 + V3) /2
    
    idxOk = abs(x-median(x))<5*std(x);
    p = polyfit(x(idxOk),y(idxOk),1);
    out(i,:) = p;
  end
end
OutTS = irf.ts_scalar(outTime,out);