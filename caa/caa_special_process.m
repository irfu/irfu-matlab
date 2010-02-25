function caa_special_process(st,dt,filename)
%CAA_SPECIAL_PROCESS: run any special processing detailed in the file
%
% caa_special_process(st,dt,filename)
% Input:
%  st,dt: start and duration (or end time) for the intervals to be
%         reprocessed, either as ISO strings or epochs.
%  filename - file with the special processing instructions
% 
% If a given interval in the file intersects the requested times, the entire
% interval is reprocessed (i.e. it's not split up by the requested times). 
%
% Example:
%   caa_special_process('2007-04-08T06:01:00Z',1, '/data/caa/l1/2007/special_processing.txt')
% Reprocesses the following entry in the special_processing.txt file:
%     2007-04-08T06:00:00Z 2007-04-09T09:00:00Z 24


if isa(st,'char'),st=iso2epoch(st); end
if isa(dt,'char'),dt=iso2epoch(dt); end
if dt>iso2epoch('2000-01-01T00:00:00Z'), dt=dt-st; end
if ~exist(filename,'file'),error(['file ' filename ' not found.']); end
et=st+dt;
old_pwd=pwd;

fid = fopen(filename);
C = textscan(fid, '%s','commentStyle', '%','delimiter','\n');
fclose(fid);
C=C{1};

for i=1:length(C)
	if ~isempty(C{i})
		entry=textscan(C{i},'%s %s %s','commentStyle', '%');
		if isISOtime(entry{1}{1}) && isISOtime(entry{2}{1})
			t1=iso2epoch(entry{1}{1});
			t2=iso2epoch(entry{2}{1});
			if t1<et && t2>st
				irf_log('proc',['Special processing interval ' C{i}])
				for cli_n=1:length(entry{3}{1})
					cli=str2double(entry{3}{1}(cli_n));
					irf_log('proc',['cli= ' num2str(cli)])
					for j=1:length(C)-i
						if isempty(C{i+j}), break, end
						entry2=textscan(C{i+j},'%s %s %s','commentStyle', '%');
						if isISOtime(entry2{1}{1}), break, end
					end
					j=j-1;
					run_special_processing(t1,t2,cli,C(i+1:i+j))
				end
			end
		end 
	end
end
cd(old_pwd);
end

function run_special_processing(st__,et__,cli,cmd__)
    old_pwd__=pwd;
	dirs__=caa_get_subdirs(epoch2iso(st__), et__-st__, cli);
	for i__=1:length(dirs__)
		cd(dirs__{i__});
		[st_int__,dt_int__] = caa_read_interval();
		st_int__=iso2epoch(st_int__);
		if (st_int__<et__) && (st_int__+dt_int__ > st__)
			irf_log('proc',['dir: ' pwd])
			for j__=1:length(cmd__)
				irf_log('proc',['cmd: ' cmd__{j__}])
				eval(cmd__{j__});
			end
		end
	end
	cd(old_pwd__);
end

function out=isISOtime(s)
  out=0;
  if length(s) < 19, return, end
  if length(s) > 30, return, end
  a = sscanf(s,'%4d-%2d-%2dT%2d:%2d:%fZ');
  if length(a) ~= 6, return, end
  if a(1) < 1900 || a(1) > 2100, return, end
  if a(2) < 1 || a(2) > 12, return, end
  if a(3) < 1 || a(3) > 31, return, end
  if a(4) < 0 || a(4) > 24, return, end
  if a(5) < 0 || a(5) > 60, return, end
  if a(6) < 0 || a(6) > 60, return, end
  out=1;
end


