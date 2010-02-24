function caa_export_month(year,month,startday,stopday)
% CAA_EXPORT_MONTH: Export all caa data for the specified month.
%
% caa_export_month(year,month,startday,stopday)
% Exports the following cef files to /data/caa/cef/YEAR:
%   L1: P1, P2, P3, P4, P12, P34, P32
%   L2: P, E
%   L3: P, E, DER
% with exceptions for probe failures.
% If stopday is not specified, takes the remainder of the month.
% If startday and stopday are not specified, defaults to the whole month.
  if nargin < 4, stopday=eomday(year,month); end
  if nargin < 3, startday=1; end

  old_pwd=pwd;
  cd(sprintf('/data/caa/cef/%4.4i',year))
  levels= [1 1 1 1 1 1 1 2 3 3 2 3];
  datatypes={'P1' 'P2' 'P3' 'P4' 'P12' 'P34' 'P32' 'P' 'P' 'DER' 'E' 'E'};
  switch year
	  case 2006
		  exceptions={{1 'P1'} {3 'P1'}  {2 'P3'} {1 'P12'} {3 'P12'} {2 'P32'} {4 'P32'}};
	  case 2007
		  if month <= 5     % probe failure on C2
			  exceptions={{1 'P1'} {2 'P3'} {3 'P1'} {1 'P12'} {2 'P32'} {3 'P12'} {4 'P32'}};
		  elseif month < 11 % software patch not yet uploaded on C2
			  exceptions={{1 'P1'} {2 'P1'} {2 'P3'} {3 'P1'} {1 'P12'} {2 'P32'} {2 'P12'} {3 'P12'} {4 'P32'}};
		  else
			  exceptions={{1 'P1'} {2 'P1'} {2 'P3'} {3 'P1'} {1 'P12'} {2 'P12'} {3 'P12'} {4 'P32'}};
		  end
	  case 2008
		  exceptions={{1 'P1'} {2 'P1'} {2 'P3'} {3 'P1'} {1 'P12'} {2 'P12'} {3 'P12'} {4 'P32'}};
	  otherwise
		  disp('Not sure which files to except.')
		  disp('To continue, enter return at the prompt.')
		  keyboard
  end
  logfile = sprintf('/data/caa/cef/%4.4i/export_log_%4.4i%2.2i.txt',year,year,month);
  irf_log('log_out',logfile);
  failfile = sprintf('/data/caa/cef/%4.4i/export_failed_%4.4i%2.2i.txt',year,year,month);
  failed=0;
  
  for day=startday:stopday
      datetxt= sprintf('%4.4i%2.2i%2.2i',year,month,day);
      disp(['--- ' datetxt ' --- '])
      start_time=toepoch([year,month,day,0,0,0]);
      for i=1:length(levels)
          level=levels(i);
          datatype=datatypes{i};
          for sat=1:4
              excepted=0;
              for j=1:length(exceptions)
                  if (exceptions{j}{1} == sat) && (strcmp(exceptions{j}{2},datatype)==1)
                      excepted=1;
                      continue
                  end
              end
              if excepted==1, continue, end
              ret=caa_export_new(level, datatype, sat, 3,'00','Null', start_time, 3600*24);
              if ret ~=0
                  fprintf(1,'caa_export_new failed for %s for lev=%1.1i. C%1.1i %s\n',datatype,level,sat,datetxt);
				  failed=1;
				  fid = fopen(failfile,'a');
				  if fid > 0
					  fprintf(fid,'caa_export_new failed: %s level %1.1i C%1.1i %s\n',datatype,level,sat,datetxt);
					  fclose(fid);
				  else
					  error(['Cannot open file' failfile])
				  end
              end
          end
      end
  end
  cd(old_pwd);
  if failed, disp(['At least one file failed to export. See ' failfile]), end
end