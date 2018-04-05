function y=irf_gradient(x,time_step)
%IRF_INTEGRATE  differenciate time series
%
% y=irf_gradient(x,time_step)
%   differenciate time series. time steps that are larger
%   than 3times the time step are assumed to be data gaps.
%
%   x - time series  to integrate
%   time_step - optional, all time_steps larger than 3*time_step are
%   assumed data gaps, default is that time_step is the second smallest 
%   value of all time_steps of the time series
%

dt=[0 ; diff(x(:,1))];
if nargin < 2 % estimate time step
    time_steps=diff(x(:,1));
    [~,ind_min]=min(time_steps);
    time_steps(ind_min)=[]; % remove the smallest time step in case some problems
    time_step=min(time_steps);
end

if size(x,2)==1, data_columns=1; else, data_columns=2:size(x,2);end
    
data_gaps=find(dt>3*time_step);
dt(data_gaps)=0;
y=x;
if data_gaps
    y(1:data_gaps(1),data_columns)=gradient(x(1:data_gaps(1),data_columns)',time_step)';
    for j=1:(length(data_gaps)-1)
        y((data_gaps(j)+1):data_gaps(j+1),data_columns)=gradient(x((data_gaps(j)+1):data_gaps(j+1),data_columns)',time_step)';
    end
    y(data_gaps(end):end,data_columns)=gradient(x(data_gaps(end):end,data_columns)',time_step)';
else
    y(:,data_columns)=gradient(x(:,data_columns)',time_step)';
end
