function caa_export_month(year,month,startday,stopday,sats,varargin)
% CAA_EXPORT_MONTH: Export all caa data for the specified month.
%
% caa_export_month(year,month,[startday],[stopday],[sats])
% Exports the following cef files to /data/caa/cef/YEAR:
%   L1: P1, P2, P3, P4, P12, P34, P32
%   L2: P, E
%   L3: P, E, DER
% with exceptions for probe failures.
% If stopday is not specified or is zero, takes the remainder of the month.
% If startday and stopday are not specified, defaults to the whole month.
% If sats is not specified, defaults to 1:4.
% Option: 'P_L23_only': exports only L2 and L3 P data.
  if nargin < 5, sats=1:4; end
  if nargin < 4, stopday=eomday(year,month); end
  if nargin < 3, startday=1; end
  if stopday==0, stopday=eomday(year,month); end 
  P_L23_only = 0;
  for i=1:length(varargin)
	switch(varargin{i})
	case 'P_L23_only'
		P_L23_only = 1;
	otherwise
		irf_log('fcal',['Option ''' varargin{i} '''not recognized'])
	end
  end

  old_pwd=pwd;
  cd(sprintf('/data/caa/cef/%4.4i',year))
  if P_L23_only
	  datatypes={ 'P' 'P' };
	  levels= [2 3];
  else
	  datatypes={'P1' 'P2' 'P3' 'P4' 'P12' 'P34' 'P32' 'P' 'P' 'DER' 'E' 'E' 'HK' 'SFIT'};
	  levels=    [1     1   1    1    1     1      1    2   3    3    2   3   2    3];
  end

  switch year
	  case 2001
		  if month <= 7   % probe 3 failure on C2
			  exceptions={ {1 'P32'} {2 'P32'} {3 'P32'} {4 'P32'}};
		  else
			  exceptions={ {2 'P3'} {1 'P32'} {2 'P32'} {3 'P32'} {4 'P32'}};
		  end
	  case 2002
		  if month <= 7   % probe 1 failure on C1
			  exceptions={ {1 'P1'} {2 'P3'} {1 'P12'} {1 'P32'} {2 'P32'} {3 'P32'} {4 'P32'}};
		  else            % probe 1 failure on C3
			  exceptions={ {1 'P1'} {2 'P3'} {3 'P1'} {1 'P12'} {3 'P12'} {1 'P32'} {2 'P32'} {3 'P32'} {4 'P32'}};
		  end
	  case 2003 % started making P32 data
		  exceptions={{1 'P1'} {3 'P1'}  {2 'P3'} {1 'P12'} {3 'P12'} {2 'P32'} {4 'P32'}};
	  case 2004
		  exceptions={{1 'P1'} {3 'P1'}  {2 'P3'} {1 'P12'} {3 'P12'} {2 'P32'} {4 'P32'}};
	  case 2005
		  exceptions={{1 'P1'} {3 'P1'}  {2 'P3'} {1 'P12'} {3 'P12'} {2 'P32'} {4 'P32'}};
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
	  case 2009
		  if month <= 10     % probe 4 failure on C1
			  exceptions={{1 'P1'} {2 'P1'} {2 'P3'} {3 'P1'} {1 'P12'} {2 'P12'} {3 'P12'} {4 'P32'}};
		  else
			  exceptions={{1 'P1'} {1 'P4'} {2 'P1'} {2 'P3'} {3 'P1'} {1 'P12'} {2 'P12'} {3 'P12'} {4 'P32'}};
		  end
	  case 2010
			  exceptions={{1 'P1'} {1 'P4'} {2 'P1'} {2 'P3'} {3 'P1'} {1 'P12'} {2 'P12'} {3 'P12'} {4 'P32'}};
	  case 2011
		  if month <= 5     % probe 3 failure on C3
			  exceptions={{1 'P1'} {1 'P4'} {2 'P1'} {2 'P3'} {3 'P1'} {1 'P12'} {2 'P12'} {3 'P12'} {4 'P32'}};
		  else
			  exceptions={{1 'P1'} {1 'P4'} {2 'P1'} {2 'P3'} {3 'P1'} {3 'P3'} {1 'P12'} {2 'P12'} {3 'P12'} {3 'P32'} {3 'P34'} {4 'P32'}};
		  end
	  otherwise
		  if year > 2011
			  disp('WARNING: uncertain which files to except.')
			  exceptions={{1 'P1'} {1 'P4'} {2 'P1'} {2 'P3'} {3 'P1'} {3 'P3'} {1 'P12'} {2 'P12'} {3 'P12'} {3 'P32'} {3 'P34'} {4 'P32'}};
		  else
			  error('Year out of range.')
		  end
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
          for sat=sats
              excepted=0;
              for j=1:length(exceptions)
                  if (exceptions{j}{1} == sat) && (strcmp(exceptions{j}{2},datatype)==1)
                      excepted=1;
                      continue
                  end
              end
              if excepted==1, continue, end
              ret=caa_export_cef(level, datatype, sat, 3,'00','Null', start_time, 3600*24);
              if ret ~=0
                  fprintf(1,'caa_export_cef failed for %s for lev=%1.1i. C%1.1i %s\n',datatype,level,sat,datetxt);
				  failed=1;
				  fid = fopen(failfile,'a');
				  if fid > 0
					  fprintf(fid,'caa_export_cef failed: %s level %1.1i C%1.1i %s\n',datatype,level,sat,datetxt);
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