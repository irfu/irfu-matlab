function caa_export_batch_l1(cl_id,out_path,lev_list)
%CAA_EXPORT_BATCH_L1 run CAA_EXPORT in a batch script
%
% caa_export_batch_l1(cl_id,[out_path],[lev_list])
%

% Copyright 2005 Yuri Khotyaintsev

if nargin<2, out_path = '.'; lev_list = 1:3;
elseif nargin==2
  if isnumeric(out_path), lev_list = out_path; out_path = '.';
  else, lev_list = 1:3;
  end
end

QUAL=3;

l1p = {'P1','P2','P3','P4'};
l1e = {'P12','P32','P34'};

if ~(exist('./mPR.mat','file') || exist('./mER.mat','file'))
  irf_log('load','no data'), return
end

v = caa_read_version;
if isempty(v), v = 0; end
v = v + 1;
v_s = num2str(v);
if v<10, v_s = ['0'  v_s]; end

fid = fopen('./.version','a');
if fid < 0, irf_log('save','problem updating version'),return, end
fprintf(fid,'%s %s\n', v_s, epoch2iso(date2epoch(now),1));
fclose(fid);

sp = pwd;
cd(out_path)

ii = find(lev_list==1);
if ~isempty(ii)
  % L1
  if exist([sp '/mPR.mat'],'file')
    for k=1:length(l1p)
      disp(['Level 1 : ' l1p{k}])
      caa_export(1,l1p{k},cl_id,QUAL,v_s,sp)
    end
  end
  if exist([sp '/mER.mat'],'file')
    for k=1:length(l1e)
      disp(['Level 1 : ' l1e{k}])
      caa_export(1,l1e{k},cl_id,QUAL,v_s,sp)
    end
  end
  lev_list(ii) = [];
end
if isempty(lev_list), cd(sp), return, end

ii = find(lev_list==2);
if ~isempty(ii)
  lev = 2;
  % EF
  if exist([sp '/mEDSIf.mat'],'file')
    disp(['Level ' num2str(lev) ' : EF'])
    caa_export(lev,'EF',cl_id,QUAL,v_s,sp)
    disp(['Level ' num2str(lev) ' : E'])
    caa_export(lev,'E',cl_id,QUAL,v_s,sp)
  end
  if exist([sp '/mP.mat'],'file')
    disp(['Level ' num2str(lev) ' : P'])
    caa_export(lev,'P',cl_id,QUAL,v_s,sp)
  end
  lev_list(ii) = [];
end

if isempty(lev_list), cd(sp), return, end

ii = find(lev_list==3);
if ~isempty(ii)
  lev = 3;
  % E
  if exist([sp '/mEDSI.mat'],'file')
    disp(['Level ' num2str(lev) ' : E'])
    caa_export(lev,'E',cl_id,QUAL,v_s,sp)
  end
  % P
  if exist([sp '/mP.mat'],'file')
    disp(['Level ' num2str(lev) ' : P'])
    caa_export(lev,'P',cl_id,QUAL,v_s,sp)
  end
end

cd(sp)
