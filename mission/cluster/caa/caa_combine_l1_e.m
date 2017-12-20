function out = caa_combine_l1_e(cl_id)
%CAA_COMBINE_L1_E  combined data from diff measureemnts into one file
% 
%  OUT = CAA_COMBINE_L1_E(CL_ID)
%
%  output: The last column contains flagP23 

out = [];

nOk = 0;
[ok,p12] = c_load('wE?p12',cl_id); if ok, nOk = nOk+1; end
[ok,p32] = c_load('wE?p32',cl_id); if ok, nOk = nOk+1; end
[ok,p34] = c_load('wE?p34',cl_id); if ok, nOk = nOk+1; end
if ~nOk, irf_log('dsrc','no data'), return, end
vars = {p12,p32,p34};

time = []; 
for v=vars
  vv = v{:};
  if ~isempty(vv), time = [time; vv(:,1)]; end %#ok<AGROW>
end; timemem=time;
time = sort(unique(time));

%%
idx = cell(length(vars),1);
if nOk > 0
  for iv=1:length(vars)
    vv = vars{iv};
    if isempty(vv), continue, end
    [~,idxTmp,~] = intersect(time,vv(:,1));
    idx{iv} = idxTmp;
  end 
end
%%
out = NaN(length(time),3);
out(:,4) = 0;
if nOk ==1
  if ~isempty(p12)
    out(:,1:2) = p12(idx{1},1:2);
  elseif ~isempty(p32)
    out(:,1:2) = p32(idx{2},1:2);
    out(:,4) = 1;
  else
    out(:,1) = p34(idx{3},1);
    out(:,3) = p34(idx{3},2);
  end
else
  out(:,1) = time;
  if ~isempty(idx{1}), out(idx{1},2) = p12(:,2); end
  if ~isempty(idx{2}), out(idx{2},2) = p32(:,2); out(idx{2},4) = 1; end
  if ~isempty(idx{3}), out(idx{3},3) = p34(:,2); end
end

end
