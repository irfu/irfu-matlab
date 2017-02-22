function out = caa_combine_l1_p(cl_id)
%CAA_COMBINE_L1_P  combined data from 4 probes into one file
% 
%  OUT = CAA_COMBINE_L1_P(CL_ID)
%

out = [];

nOk = 0;
vars = cell(1,4);
for i=1:4
  [ok,vars{i}] = c_load(irf_ssub('P10Hz?p!',cl_id,i)); if ok, nOk = nOk+1; end
end
if ~nOk, irf_log('dsrc','no data'), return, end

%%
time = []; 
for v=vars
  vv = v{:};
  if ~isempty(vv), time = [time; vv(:,1)]; end %#ok<AGROW>
end
time = sort(unique(time));

%%
idx = cell(1,length(vars));
if nOk >= 1
  for iv=1:length(vars)
    vv = vars{iv};
    if isempty(vv), continue, end
    [~,idxTmp,~] = intersect(time,vv(:,1));
    idx{iv} = idxTmp;
  end
end
%%
out = NaN(length(time),5);
out(:,1) = time;
for i=1:4
  if ~isempty(vars{i})
    out(idx{i},i+1) = vars{i}(:,2);
  end
end

end
