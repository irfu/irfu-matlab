function t_out = irf_epoch(t_in,flag)
%IRF_EPOCH  Convert time between different formats
%
%   time=IRF_EPOCH([year month date hour min sec])
%            convert [year month date hour min sec] column vector to epoch
%   time_vector=IRF_EPOCH(epoch,'vector')
%           convert epoch to [year month date hour min sec] column vector
%   time_iso=IRF_EPOCH(epoch,'iso')
%           convert to ISO format time, argument can be time in epoch
%   time_iso=IRF_EPOCH(epoch,'isoshort')
%           convert to ISO short format time, argument can be time in epoch
%
%  epoch - seconds since the epoch 1 Jan 1970.
%   The seconds since epoch time format is the time specification
%   used by the ISDAT system. To convert into Matlabs standard time
%   format use DATENUM.
%
% $Id$

if nargin==1,
    flag='toepoch';
end

switch lower(flag)
    case 'toepoch'
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
    case 'vector'
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
    case {'iso','isoshort'}
        if strcmp(flag,'isoshort'), fmt=1;end
        if length(t_in)<5
            % We need to do all this because DATESTR rounds seconds
            d = irf_epoch(t_in,'vector');
            
            for j=2:5, s1(j-1) = {add_zero(d(:,j),num2str(d(:,j),'%d'))}; end
            
            % Take care about seconds separately
            s2 = add_zero(d(:,6),num2str(d(:,6),'%6f'));
            if strcmp(flag,'isoshort'), s2 = s2(:,1:6); end
            
            sZ = s2(:,1); sZ(:) = 'Z';
            sT = sZ; sT(:) = 'T';
            sdash = sZ; sdash(:) = '-';
            scol = sZ; scol(:) = ':';
            
            t_out = [num2str(d(:,1)) sdash s1{1} sdash s1{2} sT s1{3} scol s1{4} scol s2 sZ];
        else
            % This approach is faster for data with many samples per minute, as we run
            % from epoch only once per minute
            out = char(zeros(length(t_in),27-fmt*3));
            out(:,[5 8]) = '-';
            out(:,11) = 'T';
            out(:,[14 17]) = ':';
            out(:,end) = 'Z';
            
            tss = irf_epoch(t_in(1),'vector');
            tee = irf_epoch(t_in(end),'vector');
            
            mins = irf_epoch([tss(1:5) 0]):60:irf_epoch([tee(1:5) 0]);
            d = irf_epoch(mins,'vector');
            
            s1 = {'', '', '', '',''}; sl=[0 4; 5 2; 8 2; 11 2; 14 2];
            j_start = 0;
            
            for j=1:5
                if d(1,j)==d(end,j)
                    ss = add_zero(d(1,j),num2str(d(1,j),'%d'));
                    for jj=1:sl(j,2), out(:,sl(j,1)+jj) = ss(jj); end
                else
                    j_start = j;
                    for jj=j:5, s1(jj) = {add_zero(d(:,jj),num2str(d(:,jj),'%d'))}; end
                    break
                end
            end            
            for j=1:length(mins)
                if j==length(mins), ii = find(t_in>=mins(j));
                else ii = find(t>=mins(j) & t<mins(j+1));
                end;
                if isempty(ii), continue, end
                if j_start
                    for kk=j_start:5
                        for jj=1:sl(kk,2), out(ii,sl(kk,1)+jj) = s1{kk}(j,jj); end
                    end
                end
                s2 = add_zero(t(ii)-mins(j),num2str(t(ii)-mins(j),'%6f'));
                if strcmp(flag,'isoshort'), s2 = s2(:,1:6); end
                t_out(ii,18:end-1) = s2;
            end
        end  
    otherwise
        
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
