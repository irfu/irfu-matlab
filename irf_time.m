function t_out = irf_time(t_in,flag)
%IRF_TIME  Convert time between different formats
%
%   t_out=IRF_TIME(t_in,'in2out');
%   
%   Input:
%       t_in: column vector of input time in format 'in'
%   Output:
%       t_out: column vector of output time in format 'out'
%
%   Formats 'in' and 'out' can be (default 'in' is 'epoch'):%
%      epoch: seconds since the epoch 1 Jan 1970, default, used by the ISDAT system. 
%     vector: [year month date hour min sec] in each row
%        iso: string in ISO format
%   isoshort: string in shorter ISO format 
%   yyyymmdd: string 
% yyyymmddhh: string 
%yyyymmddhhmm: string 
%       date: MATLAB datenum format
%    datenum: -=-
%        doy: [year, doy]
%         et: J2000 Ephemeris Time, seconds past  January 1, 2000, 11:58:55.816 (UTC)
%   cdfepoch: miliseconds since 1-Jan-0000
% cdfepoch16: [seconds since 1-Jan-0000, picoseconds within the second]
%
%  t_out=IRF_TIME(t_in,'out') equivalent to t_out=IRF_TIME(t_in,'epoch2out');
%  t_out=IRF_TIME(t_in) equivalent to t_out=IRF_TIME(t_in,'vector2epoch');
%
%  Example: t_out=irf_time('2011-09-13T01:23:19.000000Z','iso2epoch');
%
%  There are also commands to convert time intervals
%   time_int=IRF_TIME(tint,'tint2iso')
%           convert time interval to iso (first column epoch start time, 2nd epoch end time)  
%
% $Id$

% 'tint2iso'

persistent tlastcall strlastcall
if nargin==0, % return string with current time (second precision)
    % datestr is slow function therefore if second has not passed use the old
    % value of string as output without calling datestr function. Increases speed!
    if isempty(tlastcall),
        tlastcall=0; % initialize
    end
    if 24*3600*(now-tlastcall)>1,
        tlastcall=now;
        strlastcall=irf_time(now,'date2yyyy-mm-dd hh:mm:ss');
    end
    t_out=strlastcall;
    return
elseif nargin==1,
    flag='vector2epoch';
end

flag_tint=strfind(flag,'tint'); % check if we work with time interval (special case)
if isempty(flag_tint),          % default transformation
    flag_2=strfind(flag,'2');   % see if flag has number 2 in it 
    if isempty(flag_2),         % if no '2' convert from epoch
        format_in='epoch';   
        format_out=flag;
        flag=[format_in '2' format_out];
    else
        format_in=flag(1:flag_2-1);
        format_out=flag(flag_2+1:end);
    end
    if strcmp(format_in,format_out) % if in and out equal return
        t_out=t_in;
        return;
    elseif ~strcmp(format_in,'epoch') && ~strcmp(format_out,'epoch')
        % if there is no epoch in the format then 
        % first convert from 'in' to 'epoch' 
        % and then from 'epoch' to 'out'
        t_temp=irf_time(t_in,[format_in '2epoch']);
        t_out=irf_time(t_temp,['epoch2' format_out]);
        return
    end
end


%
% flag should include 'tint' or 'epoch'
% no other flags are allowed below
%
switch lower(flag)
    case 'vector2epoch'
        % TOEPOCH - Convert a [YYYY MM DD hh mm ss] time specification
        % to seconds since 1970.
        x=t_in;
        [m,n]=size(x);
        if n~=2 && n~=3 && n~=6,
            if m==2 || m==3 || n==6,
                x=x'; n=size(x,2);
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
    case 'epoch2vector'
        t_in=double(t_in(:));
        t = datevec(irf_time(fix(t_in),'epoch2date'));
        % THE HACK BELOW IS COMMENTED OUT! IF THERE ARE 0.01s PRECISION
        % PROBLEMS, THEY NEE TO BE UNDERSTOOD!!!!!
        % The following lines are needed to work aroung the problem with numerical
        % accuracy in conversion from isdat epoch to matlab date.
        % We give whole seconds to datevec(epoch2date)) and expect whole seconds in
        % return. Bu we get something different, and that is why we use round as we
        % expect the error to be on a level of .01 sec.
        t(:,6) = round(t(:,6));
        ii = find(t(:,6)==60);
        if ~isempty(ii)
            t(ii,6) = 0;
            t_tmp = datevec(irf_time(fix(t_in(ii)+1),'epoch2date'));
            t(ii,1:5) = t_tmp(:,1:5);
        end
        % Correct fractions of second. This actually preserves
        % accuracy ~1e-6 sec for year 2004.
        t(:,6) = t(:,6) + t_in - fix(t_in);
        t_out = t;
    case {'epoch2iso','epoch2isoshort'}
        d = irf_time(t_in,'vector');
        if strcmp(flag,'epoch2isoshort')
			fmt='%04d-%02d-%02dT%02d:%02d:%06.3fZ';
            t_out=num2str(d,fmt);
			ii=find(t_out(:,18)=='6'); % in case there has been rounding leading to 60.000 seconds
			if any(ii),
				t_out(ii,:)=num2str(irf_time(t_in(ii)+0.0005,'vector'),fmt);
			end
		else
			fmt='%04d-%02d-%02dT%02d:%02d:%09.6fZ';
            t_out=num2str(d,fmt);
			ii=find(t_out(:,18)=='6'); % in case there has been rounding leading to 60.000 seconds
			if any(ii),
				t_out(ii,:)=num2str(irf_time(t_in(ii)+0.0000005,'vector'),fmt);
			end
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
        
    case {'epoch2date','epoch2datenum'} % matlab date
        % 719529 is the number of days from 0-Jan-0000 to 1-Jan-1970
        t_out = double(719529 + double(double(t_in(:))/double(24 * 3600)));
    case {'date2epoch','datenum2epoch'} 
        t_out = double(t_in(:) - 719529)*double(24 * 3600);
        
    case 'epoch2yyyymmdd'
        t=irf_time(t_in,'epoch2vector');
        t_out=num2str(t(:,1:3),'%04d%02d%02d');
        
    case 'epoch2yyyymmddhh'
        t=irf_time(t_in,'epoch2vector');
        t_out=num2str(t(:,1:3),'%04d%02d%02d%02d');
        
    case 'epoch2yyyymmddhhmm'
        t=irf_time(t_in,'epoch2vector');
        t_out=num2str(t(:,1:5),'%04d%02d%02d%02d%02d');
                  
    case 'epoch2yyyy-mm-dd hh:mm:ss'
        d=irf_time(fix(t_in),'epoch2vector');
        t_out=num2str(d,'%04d-%02d-%02d %02d:%02d:%02.0f');
                  
    case 'epoch2doy'
          t_first_january_vector=irf_time(t_in,'vector');
          t_first_january_vector(:,2:end)=1;
          t_out=[t_first_january_vector(:,1) floor(irf_time(t_in,'datenum'))-floor(irf_time(t_first_january_vector,'vector2datenum'))+1];
          
    case 'doy2epoch'
        t_out=irf_time([t_in(:,1) t_in(:,1).*0+1 t_in(:,1).*0+1 ...
            t_in(:,2).*24-12 t_in(:,1).*0 t_in(:,1).*0]);

    case 'epoch2et'
        t_out=t_in-irf_time([2000 01 01 11 58 55.816]);

    case 'et2epoch'
        t_out=t_in+irf_time([2000 01 01 11 58 55.816]);

    case 'cdfepoch2epoch'
        t_out = (t_in-62167219200000)/1000;

    case 'epoch2cdfepoch'
        t_out = t_in*1000+62167219200000;

    case 'cdfepoch162epoch'
        t_out = (t_in(:,1)-62167219200)+t_in(:,2)*1e-12;

    case 'epoch2cdfepoch16'
        t_out = [fix(t_in)+62167219200 , mod(t_in,1)*1e12];

%
% Time interval conversions
%

    case 'tint2iso'
        t1iso=irf_time(t_in(:,1),'epoch2iso');
        t2iso=irf_time(t_in(:,2),'epoch2iso');
        t_out=[t1iso repmat('/',size(t1iso,1),1) t2iso];

    case 'iso2tint'
        % assume column array where each row is interval in iso format
        ii=strfind(t_in(1,:),'/');
        t1=irf_time(t_in(:,1:ii-1),'iso2epoch');
        t2=irf_time(t_in(ii+1:end),'iso2epoch');
        t_out=[t1 t2];

    case 'tint2isoshort'
        t1iso=irf_time(t_in(:,1),'epoch2isoshort');
        t2iso=irf_time(t_in(:,2),'epoch2isoshort');
        t_out=[t1iso repmat('/',size(t1iso,1),1) t2iso];

    otherwise
        disp(['!!! irf_time: unknown flag ''' lower(flag) ''', not converting.'])
        t_out=t_in;
end
