function caa_apply_tcor(iso_t, dt, cl_id, tcor_file)
% CAA_APPLY_TCOR: Applies the TCOR time corrections from tcor_file to all
%                 data in the interval from iso_t to iso_t+dt.
%
% caa_apply_tcor(iso_t, dt, cl_id)
%  iso_t: start time (ISO format)
%  dt:    length of interval (seconds)
%  cl_id: spacecraft to which  time correction should be applied.
%  tcor_file: TCOR file
%    (e.g. '/data/cluster/cal/WecEfwTcor/C4_CP_DWP_TCOR__20080701_000000_20090101_000000_V091022.cef')
%
% NOTE: the TCOR correction is usually performed by ISDAT, in which case
% this function should not be re-applied by calling this routine.


if ~exist(tcor_file,'file'),error(['file ' tcor_file  ' not found.']); end
fid = fopen(tcor_file,'r');
C = textscan(fid, '%s',3,'delimiter','/n');
if str2double(C{1}{1}(3)) ~= cl_id
	fclose(fid);
	error('Mismatch between cl_id and CEF file.');
end
disp('Reading TCOR file.')
C = textscan(fid, '%s %n %n','delimiter',',','CommentStyle','!');
fclose(fid);
tcor_time=iso2epoch(char(C{1})');
tcor_offset= C{2};
tcor_diff= C{3};

disp('Building directory list.')
dir_list = caa_get_subdirs(iso_t, dt, cl_id);
old_pwd=pwd;
var='';
for dir_i=1:length(dir_list)
	cd(dir_list{dir_i});
	disp(dir_list{dir_i});
	files=dir('*.mat');
	for files_i=1:length(files)
		load(files(files_i).name);
		vars=whos('-file' , files(files_i).name);
		for vars_i=1:length(vars)
			if strcmp(vars(vars_i).class,'double') && vars(vars_i).size(1) > 1 && vars(vars_i).size(2) > 1
				eval(['var=' vars(vars_i).name ';']);
				if var(1,1) > 946684800 % iso2epoch('2000-01-01T00:00:00Z')
					sz=size(var);
					if sz(2) == 2 && var(1,2) > 946684800,is_interval=1;else, is_interval=0;end
					tcor_i=caa_locate_bisect(tcor_time,var(1,1),1);
					i1=1;
					i2=0;
					valid_indx=find(isfinite(var(:,1)));
					while i2 < length(valid_indx)
						if tcor_i < length(tcor_time)
							i2=caa_locate_bisect(var(valid_indx,1),tcor_time(tcor_i+1),1);
						else
							i2=length(valid_indx);
						end
						if tcor_offset(tcor_i) > -1e29 && tcor_diff(tcor_i) > -1e29 
							if is_interval
								var(valid_indx(i1:i2),:)=var(valid_indx(i1:i2),:)+(tcor_offset(tcor_i)+tcor_diff(tcor_i))*1e-6;
							else
								var(valid_indx(i1:i2),1)=var(valid_indx(i1:i2),1)+(tcor_offset(tcor_i)+tcor_diff(tcor_i))*1e-6;
							end
						end
						i1=i2+1;
						tcor_i=tcor_i+1;
					end
				end
				eval([vars(vars_i).name '=var;']);
			end
		end
		eval(['save(''' files(files_i).name ''' ' sprintf(',''%s'' ',vars.name) ')']);
		disp(['      save(''' files(files_i).name ''' ' sprintf(',''%s'' ',vars.name) ')']);
	end
end

cd(old_pwd);
