function caa_export_batch_lev(fname, out_path, lev_list)
%CAA_EXPORT_BATCH_LEV run CAA_EXPORT in a batch script, for given levels
%
% caa_export_batch_lev(fname,[out_path],[lev_list])

% Copyright 2008 Mikael Lundberg

narginchk(1, 3)

if nargin<2, out_path = '.'; lev_list = 2:3;
elseif nargin==2
  if isnumeric(out_path), lev_list = out_path; out_path = '.';
  else, lev_list = 2:3;
  end
end

old_pwd = pwd;

QUAL=3;
DATA_VER = '00';
BASE_DIR = '/data/caa/l1';
LOG_DIR = '~/Matlab/Cluster/devel/output/export/log';
%BASE_OUTPUT = '/data/caa/cef';
BASE_OUTPUT = '~/Matlab/Cluster/devel/output/export';

time_length = 24*3600;   % Operate on 24-hour intervals.

dirs = textread(fname, '%s');
if isempty(dirs), disp('No directories given'), return, end

%ind = regexp(fname, '\.', 'once');
ind = regexp(fname, '200[1-7]\d+');
if ~isempty(ind), job_name = fname(ind(end):end);
elseif ~isempty(find(fname == '/'))
  job_name = fname(find(fname == '/', 1, 'last')+1:end);
end
disp(['----- Running job name:   ' job_name])
out_path = [BASE_OUTPUT '/' job_name];

logfile = [LOG_DIR '/' job_name '.log'];
errlog = [LOG_DIR '/' job_name '_errors'];



l1p = {'P1','P2','P3','P4'};
l1e = {'P12','P32','P34'};
l2 = {'E'};
l3 = {'E', 'DER'};
%
%if ~(exist('./mPR.mat','file') | exist('./mER.mat','file'))
%	irf_log('load','no data'), return
%end
%
%v = caa_read_version;
%if isempty(v), v = 0; end
%v = v + 1;
%v_s = num2str(v);
%if v<10, v_s = ['0'  v_s]; end
%
%fid = fopen('./.version','a');
%if fid < 0, irf_log('save','problem updating version'),return, end
%fprintf(fid,'%s %s\n', v_s, epoch2iso(date2epoch(now),1));
%fclose(fid);


dirfind = @(y) ([]);
logstatus = [];
if exist(logfile, 'file')
  logfid = fopen(logfile, 'r');
  logtext = textscan(logfid, '%s%s');
  fclose(logfid);
  [logdirs, logstatus] = deal(logtext{:});
  %   [logdirs, logstatus] = textread(logfile, '%s%s');
  if ~isempty(logdirs)
    dirfind = @(y) (find(cellfun(@(x) (~isempty(regexp(x, y))), logdirs)));
  end
end


v_s = DATA_VER;
all_errors = {};

if ~exist(out_path, 'dir')
  mkdir(out_path);
end
cd(out_path)

for kk = 1:length(dirs)
  numErrors = 0;
  d = dirs{kk};
  if regexp(d, '^%'), continue, end
  if length(d) > length(BASE_DIR) & strcmp(BASE_DIR, d(1:length(BASE_DIR)))
    d = d((length(BASE_DIR)+2):end);
  end
  cur_dir = [BASE_DIR '/' d];
  found = dirfind(d);
  if any(found & strcmp(logstatus(found), 'Done'))
    disp(['----- Skipping directory:   ' cur_dir])
    continue
  end
  disp(['----- Running directory:   ' cur_dir])

  ind = regexp(d, '_');   % Find divider between date and time in dirname.
  if isempty(ind), disp(['Invalid dir: ' d]), continue, end
  YYYY  = str2num(d(ind-8:ind-5));
  MM    = str2num(d(ind-4:ind-3));
  DD    = str2num(d(ind-2:ind-1));
  hh    = str2num(d(ind+1:ind+2));
  %   mm = str2num(d(ind+3:ind+4));
  if hh ~= 0, disp(['----- Skipping directory:   ' cur_dir]), continue, end
  start_time = toepoch([YYYY MM DD 00 00 00]);

  ind = regexp(d, 'C[1-4]'); % Find Cluster ID in dirname.
  if isempty(ind), disp(['Invalid dir: ' d]), continue, end
  cl_id = str2num(d(ind+1));

  sp = dir([cur_dir '/*_*']);

  %   keyboard

  ii = find(lev_list==1);
  if ~isempty(ii)
    % L1
    %   	if exist([sp '/mPR.mat'],'file')
    if any(checkFileExists_subfunc('mPR.mat', cur_dir, sp))
      for k=1:length(l1p)
        disp(['Level 1 : ' l1p{k}])
        if caa_export_cef(1,l1p{k},cl_id,QUAL,v_s,cur_dir,start_time,time_length) > 0
          fid = fopen(errlog, 'a');
          if fid < 0, error(['Export returned error for ' l1p{k} ' in ' cur_dir]), end
          fprintf(fid, '%s\n', ['Export returned error for ' l1p{k} ' in ' cur_dir]);
          fclose(fid);
        end
      end
    end
    %   	if exist([sp '/mER.mat'],'file')
    if any(checkFileExists_subfunc('mER.mat', cur_dir, sp))
      for k=1:length(l1e)
        disp(['Level 1 : ' l1e{k}])
        if caa_export_cef(1,l1e{k},cl_id,QUAL,v_s,cur_dir,start_time,time_length) > 0
          fid = fopen(errlog, 'a');
          if fid < 0, error(['Export returned error for ' l1e{k} ' in ' cur_dir]), end
          fprintf(fid, '%s\n', ['Export returned error for ' l1e{k} ' in ' cur_dir]);
          fclose(fid);
        end
      end
    end
    %   	lev_list(ii) = [];
  end
  %   if isempty(lev_list), cd(sp), return, end
  %   if isempty(lev_list), cd(old_pwd), return, end

  ii = find(lev_list==2);
  if ~isempty(ii)
    lev = 2;
    % EF
    %   	if exist([sp '/mEDSIf.mat'],'file')
    %      if any(checkFileExists_subfunc('mEDSIf.mat', cur_dir, sp))
    %         if findCellString_subfunc(l2, 'EF')
    %   		   disp(['Level ' num2str(lev) ' : EF'])
    %   		   caa_export_cef(lev,'EF',cl_id,QUAL,v_s,cur_dir,start_time,time_length)
    %         end
    if findCellString_subfunc(l2, 'E')
      disp(['Level ' num2str(lev) ' : E'])
      try
        if caa_export_cef(lev,'E',cl_id,QUAL,v_s,cur_dir,start_time,time_length) > 0
          fid = fopen(errlog, 'a');
          if fid < 0, error(['Export returned error for L2E in ' cur_dir]), end
          fprintf(fid, '%s\n', ['Export returned error for L2E in ' cur_dir]);
          fclose(fid);
          numErrors = numErrors + 1;
        end
      catch
        logfid = fopen(logfile, 'a');
        if logfid < 0, warning(['Cannot open log file: ' logfile '. Proceeding without logging.']), end
        fprintf(logfid, '%s\tError_L2E\n', cur_dir);
        fclose(logfid);
        numErrors = numErrors + 1;
      end
      %   		end
    end
    %   	if exist([sp '/mP.mat'],'file')
    if any(checkFileExists_subfunc('mP.mat', cur_dir, sp))
      if findCellString_subfunc(l2, 'P')
        disp(['Level ' num2str(lev) ' : P'])
        if caa_export_cef(lev,'P',cl_id,QUAL,v_s,cur_dir,start_time,time_length) > 0
          fid = fopen(errlog, 'a');
          if fid < 0, error(['Export returned error for L2P in ' cur_dir]), end
          fprintf(fid, '%s\n', ['Export returned error for L2P in ' cur_dir]);
          fclose(fid);
        end
      end
    end
    %   	lev_list(ii) = [];
  end
  %   if isempty(lev_list), cd(sp), return, end
  %   if isempty(lev_list), cd(old_pwd), return, end

  ii = find(lev_list==3);
  if ~isempty(ii)
    lev = 3;
    % E
    %   	if exist([sp '/mEDSI.mat'],'file')
    %     if any(checkFileExists_subfunc('mEDSI.mat', cur_dir, sp))
    if findCellString_subfunc(l3, 'E')
      disp(['Level ' num2str(lev) ' : E'])
      try
        if caa_export_cef(lev,'E',cl_id,QUAL,v_s,cur_dir,start_time,time_length) > 0
          fid = fopen(errlog, 'a');
          if fid < 0, error(['Export returned error for L3E in ' cur_dir]), end
          fprintf(fid, '%s\n', ['Export returned error for L3E in ' cur_dir]);
          fclose(fid);
        end
      catch
        logfid = fopen(logfile, 'a');
        if logfid < 0, warning(['Cannot open log file: ' logfile '. Proceeding without logging.']), end
        fprintf(logfid, '%s\tError_L3E\n', cur_dir);
        fclose(logfid);
        numErrors = numErrors + 1;
      end
      %   		end
      if findCellString_subfunc(l3, 'DER')
        disp(['Level ' num2str(lev) ' : DER'])
        try
          if caa_export_cef(lev,'DER',cl_id,QUAL,v_s,cur_dir,start_time,time_length) > 0
            fid = fopen(errlog, 'a');
            if fid < 0, error(['Export returned error for L3DER in ' cur_dir]), end
            fprintf(fid, '%s\n', ['Export returned error for L3DER in ' cur_dir]);
            fclose(fid);
          end
        catch
          logfid = fopen(logfile, 'a');
          if logfid < 0, warning(['Cannot open log file: ' logfile '. Proceeding without logging.']), end
          fprintf(logfid, '%s\tError_L3DER\n', cur_dir);
          fclose(logfid);
          numErrors = numErrors + 1;
        end
      end
    end
    % P
    %   	if exist([sp '/mP.mat'],'file')
    if any(checkFileExists_subfunc('mP.mat', cur_dir, sp))
      if findCellString_subfunc(l3, 'P')
        disp(['Level ' num2str(lev) ' : P'])
        if caa_export_cef(lev,'P',cl_id,QUAL,v_s,cur_dir,start_time,time_length) > 0
          fid = fopen(errlog, 'a');
          if fid < 0, error(['Export returned error for L3P in ' cur_dir]), end
          fprintf(fid, '%s\n', ['Export returned error for L3P in ' cur_dir]);
          fclose(fid);
        end
      end
    end
  end   % lev == 3

  if numErrors == 0
    logfid = fopen(logfile, 'a');
    if logfid < 0, warning(['Cannot open log file: ' logfile '. Proceeding without logging.']), end
    fprintf(logfid, '%s\tDone\n', cur_dir);
    fclose(logfid);
  end
end   % for
cd(old_pwd)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subfunction to check whether a file 'fn' exists in a
% given list of subdirs 'sp' do a directory 'cur_dir'
function found = checkFileExists_subfunc(fn, cur_dir, sp)
found = [];
for ff = 1:length(sp)
  found(ff) = exist([cur_dir '/' sp(ff).name '/' fn], 'file');
end

% Subfunction to find if a string 'theString' exists in
% a cell array of strings 'theCell'. Returns index if found.
function found = findCellString_subfunc(theCell, theString)
for k = 1:length(theCell)
  if (theCell{k} == theString), found = k; return, end
end
found = 0;