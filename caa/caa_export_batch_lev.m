function caa_export_batch_lev(fname, out_path, lev_list)
%CAA_EXPORT_BATCH_LEV run CAA_EXPORT in a batch script, for given levels
%
% caa_export_batch_lev(fname,[out_path],[lev_list])
%
% $Id$

% Copyright 2008 Mikael Lundberg

error(nargchk(1, 3, nargin))

if nargin<2, out_path = '.'; lev_list = 1:3; 
elseif nargin==2
	if isnumeric(out_path), lev_list = out_path; out_path = '.';
	else, lev_list = 1:3;
	end
end

old_pwd = pwd;

QUAL=3;
DATA_VER = '05';
BASE_DIR = '/data/caa/l1';
BASE_OUTPUT = '/data/caa/cef';
time_length = 3*3600;   % Operate on 3-hour intervals.

dirs = textread(fname, '%s');
if isempty(dirs), disp('No directories given'), return, end




l1p = {'P1','P2','P3','P4'};
l1e = {'P12','P32','P34'};
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



v_s = DATA_VER;

if ~exist(out_path, 'dir')
   mkdir(out_path);
end
cd(out_path)

for kk = 1:length(dirs)
   d = dirs{kk};
   if length(d) > length(BASE_DIR) & strcmp(BASE_DIR, d(1:length(BASE_DIR)))
      d = d((length(BASE_DIR)+2):end);
   end
   cur_dir = [BASE_DIR '/' d];
   disp(['Running directory:   ' cur_dir])
   keyboard
   ind = regexp(d, '_');   % Find divider between date and time in dirname.
   if isempty(ind), disp(['Invalid dir: ' d]), continue, end
   YYYY  = str2num(d(ind-8:ind-5));
   MM    = str2num(d(ind-4:ind-3));
   DD    = str2num(d(ind-2:ind-1));
   hh    = str2num(d(ind+1:ind+2));
%   mm = str2num(d(ind+3:ind+4));
   start_time = toepoch([YYYY MM DD hh 00 00]);
   
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
   			caa_export(1,l1p{k},cl_id,QUAL,v_s,cur_dir,start_time,time_length) 
   		end
   	end
%   	if exist([sp '/mER.mat'],'file')
      if any(checkFileExists_subfunc('mER.mat', cur_dir, sp))
   		for k=1:length(l1e)
   			disp(['Level 1 : ' l1e{k}])
   			caa_export(1,l1e{k},cl_id,QUAL,v_s,cur_dir,start_time,time_length) 
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
      if any(checkFileExists_subfunc('mEDSIf.mat', cur_dir, sp))
   		disp(['Level ' num2str(lev) ' : EF'])
   		caa_export(lev,'EF',cl_id,QUAL,v_s,cur_dir,start_time,time_length)
   		disp(['Level ' num2str(lev) ' : E'])
   		caa_export(lev,'E',cl_id,QUAL,v_s,cur_dir,start_time,time_length)
   	end
%   	if exist([sp '/mP.mat'],'file')
     if any(checkFileExists_subfunc('mP.mat', cur_dir, sp))
   		disp(['Level ' num2str(lev) ' : P'])
   		caa_export(lev,'P',cl_id,QUAL,v_s,cur_dir,start_time,time_length)
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
     if any(checkFileExists_subfunc('mEDSI.mat', cur_dir, sp))
   		disp(['Level ' num2str(lev) ' : E'])
   		caa_export(lev,'E',cl_id,QUAL,v_s,cur_dir,start_time,time_length)
   		disp(['Level ' num2str(lev) ' : DER'])
   		caa_export(lev,'DER',cl_id,QUAL,v_s,cur_dir,start_time,time_length)
   	end
   	% P
%   	if exist([sp '/mP.mat'],'file')
     if any(checkFileExists_subfunc('mP.mat', cur_dir, sp))
   		disp(['Level ' num2str(lev) ' : P'])
   		caa_export(lev,'P',cl_id,QUAL,v_s,cur_dir,start_time,time_length)
   	end
   end
end
cd(old_pwd)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function found = checkFileExists_subfunc(fn, cur_dir, sp)
   for ff = 1:length(sp)
      found(ff) = exist([cur_dir '/' sp(ff).name '/' fn], 'file');
   end