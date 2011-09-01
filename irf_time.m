function t_out = irf_time(t_in,flag)
%IRF_TIME  Convert time between different formats
%
%   epoch=IRF_TIME([year month date hour min sec])
%   epoch=IRF_TIME([year month date hour min sec],'vector2epoch')
%            convert [year month date hour min sec] column vector to epoch
%
%   time_vector=IRF_TIME(epoch,'epoch2vector')
%   time_vector=IRF_TIME(epoch,'vector')
%           convert epoch to [year month date hour min sec] column vector
%
%   time_iso=IRF_TIME(epoch,'epoch2iso')
%   time_iso=IRF_TIME(epoch,'iso')
%           convert to ISO format time, argument can be time in epoch
%
%   epoch=IRF_TIME(time_iso,'iso2epoch')
%
%   time_iso=IRF_TIME(epoch,'epoch2isoshort')
%   time_iso=IRF_TIME(epoch,'isoshort')
%           convert to ISO short format time
%
%   time_yyyymmdd=IRF_TIME(epoch,'epoch2yyyymmdd')
%   time_yyyymmdd=IRF_TIME(epoch,'yyyymmdd')
%
%   date=IRF_TIME(epoch,'epoch2date')
%   date=IRF_TIME(epoch,'date')
%           convert to Matlab date format 
%
%   epoch=IRF_TIME(date,'date2epoch')
%           convert from Matlab date to epoch 
%
%   doy=IRF_TIME(date,'date2doy')
%   date=IRF_TIME([year doy],'doy2date')
%           convert between date [yyyy mm dd hh mm ss] and doy 
%
%   time_int=IRF_TIME(tint,'tint2iso')
%           convert time interval to iso (first column epoch start time, 2nd epoch end time)  
%
%  epoch - seconds since the epoch 1 Jan 1970.
%   The seconds since epoch time format is the time specification
%   used by the ISDAT system. 
%
% $Id$

% 'tint2iso'

if nargin==1,
    flag='vector2epoch';
end

switch lower(flag)
    case 'vector2epoch'
        % TOEPOCH - Convert a [YYYY MM DD hh mm ss] time specification
        % to seconds since 1970.
        x=t_in;
        [m,n]=size(x);
        if n~=2 && n~=3 && n~=6,
            if m==2 || m==3 || n==6,
                x=x'; [~,n]=size(x);
            else
                warning('irfu:argument','irf_time:Illegal argument')
                t_out=NaN;
                return
            end
        end
        
        if n==2,
            y=x;
            x(:,[3,6])=rem(y,100);
            x(:,[2,5])=rem(floor(y/100),100);
            x(:,[1,4])=floor(y/10000);
        elseif n==3,
            x(:,6)=rem(x(:,3),100); x(:,5)=floor(x(:,3)/100);
            x(:,4)=rem(x(:,2),100); x(:,3)=floor(x(:,2)/100);
            x(:,2)=rem(x(:,1),100); x(:,1)=floor(x(:,1)/100)+1900;
        elseif n==6,
            if x(1,1)<100, x(:,1)=1900+x(:,1); end
        end
        
        years=x(:,1);hours=zeros(size(years));secs=hours;
        for year=unique(x(:,1))'
            daym=[0 31 28 31 30 31 30 31 31 30 31 30 31];
            if rem(year,4)==0, daym(3)=29; end  % works up to 2100
            days=cumsum(daym(1:12))';
            
            ind=find(years==year);
            hours(ind,1)=(days(x(ind,2))+(x(ind,3)-1))*24+x(ind,4);
            secs(ind,1)=(hours(ind)*60+x(ind,5))*60+x(ind,6);
            plus=0;
            diff_yr = year-1970;
            for i = 1:diff_yr
                if rem(1969+i,4)==0
                    plus = plus + 31622400;
                else
                    plus = plus + 31536000;
                end
            end
            secs(ind,1)=secs(ind,1)+plus;
        end
        t_out=secs;
    case {'epoch2vector','vector'}
        t = datevec(irf_time(fix(double(t_in(:))),'epoch2date'));
        % The following lines are needed to work aroung the problem with numerical
        % accuracy in conversion from isdat epoch to matlab date.
        % We give whole seconds to datevec(epoch2date)) and expect whole seconds in
        % return. Bu we get something different, and that is why we use round as we
        % expect the error to be on a level of .01 sec.
        t(:,6) = round(t(:,6));
        ii = find(t(:,6)==60);
        if ~isempty(ii)
            t(ii,6) = 0;
            t_tmp = datevec(irf_time(fix(double(t_in(ii)+1)),'epoch2date'));
            t(ii,1:5) = t_tmp(:,1:5);
        end
        % Correct fractions of second. This actually preserves
        % accuracy ~1e-6 sec for year 2004.
        t(:,6) = t(:,6) + double(t_in(:)) - fix(double(t_in(:)));
        t_out = t;
    case {'iso','isoshort','epoch2iso','epoch2isoshort'}
        d = irf_time(t_in,'vector');
        if strcmp(flag,'isoshort')|| strcmp(flag,'epoch2isoshort')
            t_out=num2str(d,'%04d-%02d-%02dT%02d:%02d:%06.3fZ');
        else
            t_out=num2str(d,'%04d-%02d-%02dT%02d:%02d:%09.6fZ');
        end
    case 'iso2epoch'
        mask = '%4d-%2d-%2dT%2d:%2d:%fZ';
        s=t_in;
        % If we have multiple rows, we need to turn the matrix
        if min(size(s))>1
            if size(s,2)==27 || size(s,2)==24, s=s'; end
            n_column = size(s,2);
        else n_column = 1;
        end
        
        a = sscanf(s,mask);
        
        N = length(a)/6;
        if N~=fix(N) || N~=n_column, disp('something is wrong with input'), end
        a = reshape(a,6,fix(N));
        a = a';
        t_out = irf_time([a(:,1) a(:,2) a(:,3) a(:,4) a(:,5) a(:,6)]);
        
    case {'date','epoch2date'} % matlab date
        % 719529 is the number of days from 0-Jan-0000 to 1-Jan-1970
        t_out = double(719529 + double(double(t_in(:))/double(24 * 3600)));
        
    case 'date2epoch' 
        t_out = double(t_in(:) - 719529)*double(24 * 3600);
        
    case {'yyyymmdd','epoch2yyyymmdd'}
        t=irf_time(t_in,'epoch2vector');
        t_out=sprintf('%04d%02d%02d',t(1),t(2),t(3));
        
    case 'date2doy'
          if nargin<1, help date2doy, return, end
          eomdays = eomday(t_in(:,1)*ones(1,12),ones(length(t_in(:,1)),1) *(1:12));
          doy=[];
          for k=1:size(t_in,1)
              t_out = [doy; sum( eomdays(k,1:t_in(k,2)-1),2)];
          end
          t_out = t_out + t_in(:,3);
          
    case 'doy2date'
        YEAR=t_in(:,1);
        DOY=t_in(:,2);
        t_out = ones(length(YEAR(:,1)),3);
        for n=1:length(YEAR(:,1))
            year = YEAR(n,:); doy  =  DOY(n,:);
            total_day = cumsum(eomday(year,1:12));
            month = find( total_day>= doy ); month = month(1);
            if month>1
                pmonth = find( total_day< doy ); day = doy-total_day(pmonth(end));
            else
                day = doy;
            end
            t_out(n,:) = [year month day];
        end
        
  case 'tint2iso'
        t1iso=irf_time(t_in(:,1),'epoch2iso');
        t2iso=irf_time(t_in(:,2),'epoch2iso');
        t_out=[t1iso repmat('/',size(t1iso,1),1) t2iso];


    otherwise
        disp('!!! irf_time: unknown flag, not converting.')
        t_out=t_in;
end
