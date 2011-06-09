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
%   time_iso=IRF_TIME(epoch,'epoch2isoshort')
%   time_iso=IRF_TIME(epoch,'isoshort')
%           convert to ISO short format time, argument can be time in epoch
%
%   time_yyyymmdd=IRF_TIME(epoch,'epoch2yyyymmdd')
%   time_yyyymmdd=IRF_TIME(epoch,'yyyymmdd')
%           convert to YYYYMMDD format
%
%   date=IRF_TIME(epoch,'epoch2date')
%   date=IRF_TIME(epoch,'date')
%           convert to Matlab date format 
%
%   epoch=IRF_TIME(date,'date2epoch')
%           convert from Matlab date to epoch 
%
%  epoch - seconds since the epoch 1 Jan 1970.
%   The seconds since epoch time format is the time specification
%   used by the ISDAT system. 
%
% $Id$

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
                warning('Illegal argument')
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
        
        years=x(:,1);
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
        t = datevec(epoch2date(fix(double(t_in(:)))));
        % The following lines are needed to work aroung the problem with numerical
        % accuracy in conversion from isdat epoch to matlab date.
        % We give whole seconds to datevec(epoch2date)) and expect whole seconds in
        % return. Bu we get something different, and that is why we use round as we
        % expect the error to be on a level of .01 sec.
        t(:,6) = round(t(:,6));
        ii = find(t(:,6)==60);
        if ~isempty(ii)
            t(ii,6) = 0;
            t_tmp = datevec(epoch2date(fix(double(t_in(ii)+1))));
            t(ii,1:5) = t_tmp(:,1:5);
        end
        % Correct fractions of second. This actually preserves
        % accuracy ~1e-6 sec for year 2004.
        t(:,6) = t(:,6) + double(t_in(:)) - fix(double(t_in(:)));
        t_out = t;
    case {'iso','isoshort','epoch2iso','epoch2isoshort'}
        if length(t_in)<5
            % We need to do all this because DATESTR rounds seconds
            d = irf_time(t_in,'vector');
            
            for j=2:5, 
                s1(j-1) = {num2str(d(:,j),'%02d')}; 
            end
            
            % Take care about seconds separately
            s2 = num2str(d(:,6),'%09.6f');
            if strcmp(flag,'isoshort') || strcmp(flag,'epoch2isoshort'), 
                s2 = num2str(d(:,6),'%06.3f');
            end
            
            sZ = s2(:,1); sZ(:) = 'Z';
            sT = sZ; sT(:) = 'T';
            sdash = sZ; sdash(:) = '-';
            scol = sZ; scol(:) = ':';
            
            t_out = [num2str(d(:,1)) sdash s1{1} sdash s1{2} sT s1{3} scol s1{4} scol s2 sZ];
        else
            % This approach is faster for data with many samples per minute, as we run
            % from epoch only once per minute
            if strcmp(flag,'isoshort') || strcmp(flag,'epoch2isoshort') , 
                fmt=1;
            else
                fmt=0;
            end
            t_out = char(zeros(length(t_in),27-fmt*3));
            t_out(:,[5 8]) = '-';
            t_out(:,11) = 'T';
            t_out(:,[14 17]) = ':';
            t_out(:,end) = 'Z';
            
            tss = irf_time(t_in(1),'vector');
            tee = irf_time(t_in(end),'vector');
            
            mins = irf_time([tss(1:5) 0]):60:irf_time([tee(1:5) 0]);
            d = irf_time(mins,'vector');
            
            s1 = {'', '', '', '',''}; sl=[0 4; 5 2; 8 2; 11 2; 14 2];
            j_start = 0;
            
            for j=1:5
                if d(1,j)==d(end,j)
                    ss = add_zero(d(1,j),num2str(d(1,j),'%d'));
                    for jj=1:sl(j,2), 
                        t_out(:,sl(j,1)+jj) = ss(jj); 
                    end
                else
                    j_start = j;
                    for jj=j:5, s1(jj) = {add_zero(d(:,jj),num2str(d(:,jj),'%d'))}; end
                    break
                end
            end            
            for j=1:length(mins)
                if j==length(mins), ii = find(t_in>=mins(j));
                else ii = find(t_in>=mins(j) & t_in<mins(j+1));
                end;
                if isempty(ii), continue, end
                if j_start
                    for kk=j_start:5
                        for jj=1:sl(kk,2), t_out(ii,sl(kk,1)+jj) = s1{kk}(j,jj); end
                    end
                end
                s2 = add_zero(t_in(ii)-mins(j),num2str(t_in(ii)-mins(j),'%6f'));
                if strcmp(flag,'isoshort')|| strcmp(flag,'epoch2isoshort')
                    s2 = s2(:,1:6); 
                end
                t_out(ii,18:end-1) = s2;
            end
        end  
    case {'date','epoch2date'} % matlab date
        t_out = double(719529 + double(double(t_in(:))/double(24 * 3600)));
    case 'date2epoch' 
        t_out = double(t_in(:) - 719529)*double(24 * 3600);
    case {'yyyymmdd','epoch2yyyymmdd'}
        t=fromepoch(t_in);
        t_out=sprintf('%04d%02d%02d',t(1),t(2),t(3));
    otherwise
        disp('!!! irf_time: unknown flag, not converting.')
        t_out=t_in;
end

% Help function to insert zeros for numbers containing only one digit
function out = add_zero(d,s)
out = s;
ii = find(d<10);
if ~isempty(ii)
    ss = s(ii,1);
    ss(:) = '0';
    if length(ii)==length(d)
        % Add to all lines
        out = [ss s];
    else
        out(ii,:) = [ss s(ii,1:end-1)];
    end
end
