function [z]=av_t_appl(x,s)
%function [z]=av_t_appl(x,s)
% x is time vector, first column time  other data
% s is what expression to apply to data part
% ex. y=av_t_appl(x,'*2/1e3');

warning('caa:cleanup',...
'Function %s is deprecated and will be removed on May 1, 2004.\nUse %s instead',...
mfilename,'irf_tappl')

z=x;
z(:,2:end)=eval(['x(:,2:end)',s]);
