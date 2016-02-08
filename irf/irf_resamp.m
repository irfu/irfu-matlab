function out = irf_resamp(x,y,varargin)
%IRF_RESAMP   Resample X to the time line of Y
%
% if sampling of X is more than two times higher than Y, we average X,
% otherwise we interpolate X.
%
% out = irf_resamp(X,Y,[METHOD],['fsample',FSAMPLE],['window',WIN],
%                      ['thresh',THRESH],['median'],['max'])
% method - method of interpolation 'spline', 'linear' etc. (default 'linear')
%          if method is given then interpolate independant of sampling
% thresh - points above STD*THRESH are disregarded for averaging
% fsample - sampling frequency of the Y signal, 1/window
% window - length of the averaging window, 1/fsample
% median - use median instead of mean when averaging
% max    - return max within each averaging window, rather than mean
%
% See also INTERP1

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

narginchk(2,8)

have_options = 0;
args = varargin; 
if nargin > 2, have_options = 1; end

% Default values that can be override by options
sfy = [];
thresh = 0;
method = '';
flag_do='check'; % if no method check if interpolate or average
median_flag=0;
max_flag=0;
%doReshape = 0;

while have_options
	l = 1;
  aa = lower(args{1});
	switch(lower(args{1}))
		case {'nearest','linear','spline','pchip','cubic','v5cubic'}
			method = args{1};
			flag_do='interpolation'; % if method is given do interpolation
		case 'method'
			if length(args)>1
				if ischar(args{2})
					method = args{2};
					l = 2;
					flag_do='interpolation'; % if method is given do interpolation
				else irf.log('critical','wrongArgType : METHOD must be numeric')
				end
			else irf.log('critical','wrongArgType : METHOD value is missing')
			end
		case {'fs','fsample'}
			if length(args)>1
				if ~isempty(sfy)
          msgS = 'FSAMPLE/WINDOW already specified';
          irf.log('critical',msgS), error(msgS)
				end
				if isnumeric(args{2})
					sfy = args{2};
					l = 2;
				else irf.log('critical','wrongArgType : FSAMPLE must be numeric')
				end
			else irf.log('critical','wrongArgType : FSAMPLE value is missing')
			end
		case {'win','window'}
			if length(args)>1
				if ~isempty(sfy)
					msgS = 'FSAMPLE/WINDOW already specified';
          irf.log('critical',msgS), error(msgS)
				end
				if isnumeric(args{2})
					sfy = 1/args{2};
					l = 2;
				else irf.log('critical','wrongArgType : WINDOW must be numeric')
				end
			else irf.log('critical','wrongArgType : WINDOW value is missing')
			end
		case {'thresh','threshold'}
            if length(args)>1
                if isnumeric(args{2})
                    thresh = args{2};
                    l = 2;
                else irf.log('critical','wrongArgType : THRESHOLD must be numeric')
                end
            else irf.log('critical','wrongArgType : THRESHOLD value is missing')
            end
        case 'median'
            median_flag=1;
        case 'max'
            max_flag=1;
		otherwise
			irf.log('warning',['Skipping parameter ''' args{1} ''''])
			args = args(2:end);
	end
	args = args(l+1:end);
	if isempty(args), break, end
end

% return in case inputs are empty
% changed numel(x) == 0 to isempty(x), since it also works on TSeries
if isempty(x) || isempty(y)
    out=[];
    irf.log('warining','Some of input is empty, returning empty ouput');
    return;
end

% Need to check that x and y are given in same or similar enough format,
if any(isa(x,'TSeries')) && ~any([isa(y,'TSeries') isa(y,'GenericTimeArray')])
  error('If X is TSeries, Y must be TSeries or GenericTimeArray.')
end
  
% Check format of input
xOrig = x; yOrig = y;

% x, old time and data
if isa(x,'TSeries')
  flag_output = 'tseries';      
  oldData = x.data;  
  oldTime = double(x.time.ttns-x.time(1).ttns);%/1e9; % transform from nanoseconds to seconds
elseif isa(x,'numeric') % old format [time data]
  flag_output = 'mat';
  if numel(x) == 1 % only one number
    oldData = x;
    oldTime = [];
  else 
    oldData = x(:,2:end);   
    oldTime = x(:,1);       
  end    
else
  msgS = 'Cannot recognize input x (old data), it is neither a TSeries nor a numeric array';
  irf.log('critical',msgS), error(msgS)
end

flag_datatype = class(oldData);
% y, new time
if isa(y,'struct'),
  if isfield(y,'t'), t=y.t; t=t(:);
  else
    msgS = 'Input is structure without time field';
    irf.log('critical',msgS), error(msgS)
  end
elseif isa(y,'TSeries')
  y = y.time;
  type_epoch = class(y);
  t = double(y.time.ttns-x.time(1).ttns);%/1e9; % transform from nanoseconds to seconds
elseif isa(y,'GenericTimeArray')
  type_epoch = class(y);
  t = double(y.ttns-x.time(1).ttns);%/1e9;  
elseif isa(y,'numeric')
  if size(y,2)==1, t = y(:);  % y is only time
  else t = y(:,1); t = t(:);  % first column of y is time
  end
else 
  msgS = 'Cannot recognize input y (new time), it is neither a TSeries nor a numeric array';
  irf.log('critical',msgS), error(msgS)
end

% Reshape data so that it becomes a 2D matrix (nTimes x nData), and reshape
% back later tot he original dimensions (except the time)
origOldData = oldData;
origDataSize = size(oldData);
oldData = reshape(oldData,[origDataSize(1) prod(origDataSize(2:end))]);

% Same timeline - no need to do anything
% old data already divided into time and data
if length(oldTime)==length(t) && all(oldTime==t)
  irf.log('notice','New and old timelines are identical - no resampling needed')
  out = x; 
  return
end

% If X (oldData) has only one point, this is a trivial case and we 
% return directly there is only one number
% the single data point is duplicated into the new timeline
if numel(oldTime) == 1,       
  out = construct_output(t,repmat(oldData,numel(t),1,1));
  return
end

ndata = length(t);
if strcmp(flag_do,'check'), % Check if interpolation or average
	if ndata>1 
		% If more than one output time check sampling frequencies
		% to decide interpolation/average
		
		% Guess samplings frequency for Y   
		if isempty(sfy)
			sfy1 = (1/(t(2) - t(1)));
			if ndata==2, sfy = sfy1;
			else
				not_found = 1; cur = 3; MAXTRY = 10;
				while (not_found && cur<=ndata && cur-3<MAXTRY)
					sfy = 1/(t(cur) - t(cur-1));
            if abs(sfy-sfy1)<sfy*0.001
                not_found = 0;
                sfy = (sfy+sfy1)/2;
                break
            end
            sfy1=sfy;
					cur = cur + 1;
				end
				if not_found
					sfy = sfy1;
					irf.log('warning',	sprintf(...
						'Cannot guess sampling frequency. Tried %d times',MAXTRY));
				end
			end
			clear sfy1
		end
		
		if length(oldTime(:,1))/(oldTime(end,1) - oldTime(1,1)) > 2*sfy
			flag_do='average';
			irf.log('warning','Using averages in irf_resamp.');
		else
			flag_do='interpolation';
		end
	else
		flag_do='interpolation';  % If one output time then do interpolation
	end
end

if strcmp(flag_do,'average')
  % Input to irf_average_mx needs to be double
  oldData = double(oldData);  
  dt2 = .5/sfy; % Half interval
  if median_flag || max_flag || (exist('irf_average_mx','file')~=3)
      if (~median_flag && ~max_flag), irf.log('warning','cannot find mex file, defaulting to Matlab code.')
      end
      newData = zeros(ndata,size(oldData,2),size(oldData,3));      
      for j=1:ndata
          ii = find(oldTime(:,1) <=  t(j) + dt2 & oldTime(:,1) >  t(j) - dt2);            
          if isempty(ii), newData(j,:) = NaN;
          else
              if thresh % Throw away points above THERESH*STD()
                  sdev = std(oldData(ii,:));
                  mm = mean(oldData(ii,:));
                  if any(~isnan(sdev))
                      for k=1:length(sdev)
                          if ~isnan(sdev(k))
                              kk = find( abs( oldData(ii,k) -mm(k) ) <= thresh*sdev(k));
                              %disp(sprintf(...
                              %	'interval(%d) : disregarding %d 0f %d points',...
                              %	j, length(ii)-length(kk),length(ii)));
                              if ~isempty(kk)
                                  if median_flag, newData(j,k) = median(oldData(ii(kk),k));
                                  elseif max_flag, newData(j,k) = max(oldData(ii(kk),k));
                                  else newData(j,k) = mean(oldData(ii(kk),k));
                                  end
                              end
                          else newdata(j,k) = NaN;
                          end
                      end
                  else newData(j,:) = NaN;
                  end
              else
                  if median_flag, newData(j,:,:) = median(oldData(ii,:,:));
                  elseif max_flag, newData(j,:,:) = max(oldData(ii,:,:));
                  else newData(j,:,:) = mean(oldData(ii,:,:));
                  end
              end
          end
      end
      out = construct_output(t,newData);
  else
      if ( (oldTime(1,1) > t(end) + dt2) || (oldTime(end,1) <= t(1) - dt2) )
          irf.log('warning','Interval mismatch - empty return')           
          out = construct_output([],[]);
      else          
        tmpData = irf_average_mx([oldTime oldData],t,double(dt2),thresh);               
        newData = tmpData(:,2:end);
      end
      out = construct_output(t,newData);
  end
elseif strcmp(flag_do,'interpolation'),
  if ~any([strcmp(flag_datatype,'double') strcmp(flag_datatype,'single')])
    % Maybe insert warning here that data is converted to double, if it is 
    % not already double
    oldData = double(oldData);
  end
  if nargin < 3 || isempty(method), method = 'linear'; end

  % If time series agree, no interpolation is necessary.
  if size(oldData,1)==size(t,1), if oldTime==t(:,1), out = construct_output(oldTime,oldData); return, end, end
  
  out = construct_output(t,interp1(oldTime,oldData,t,method,'extrap'));
end

function out = construct_output(t,newData)
% Constructing output from t/y and newData
% No real need to pass t?
  if isempty(t) || isempty(newData)
    out = [];
    return;
  end  
  
  % Reshape data back so that it regains its old dimensions  
  newData = reshape(newData,[length(t) origDataSize(2:end)]); % shape back to original dimensions
  
  switch flag_output
    case 'tseries'              
      newData = cast(newData,flag_datatype);  % Recast data into original datatype
      out = xOrig.clone(y,newData);              
    case 'mat'
      out = [t newData];
  end
end % end of construct_output()

end % end of irf_resamp()
